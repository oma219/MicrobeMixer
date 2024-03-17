
# Name: simulate.py
# Description: This code contains methods related to simulate sub-command
#              in the main script.

from stats_requests import parse_taxonomy_file
from utils import log_message, run_command
from multiprocessing import Pool
import regex
import random
import os

def simulate_main(args):
    """ take the top n most abundant genera and simulate 16S rRNA reads"""
    
    # get a list of the genera from the taxonomy
    genera_to_count, _ = parse_taxonomy_file(args.silva_tax_path, ignore_eukaryotes=False)
    genera_to_seqs = {genus: [] for genus in genera_to_count.keys()}

    # process the reference file and store each sequence with its corresponding genus
    genera_to_seqs = process_reference_file(args.silva_ref_path, genera_to_seqs)

    # process the forward and reverse primer sequences
    for_primer, rev_primer = process_primer_file(args.primer_file)

    # load the top n most abundant genera file
    top_genera_to_stats = process_biome_abundance_file(args.biome_abund_path)
    assert all([genus in genera_to_seqs for genus in top_genera_to_stats.keys()]), "not all genera in top n file are in reference file"
    
    # extract the variable regions using in-silico PCR for each top n genus
    top_genera_extracted_regions = {top_genus: [] for top_genus in top_genera_to_stats.keys()}
    top_genera_extracted_regions, empty_genera = extract_variable_regions(top_genera_extracted_regions, genera_to_seqs, for_primer, rev_primer)
    
    # write out the extracted sequences to reference files
    top_genera_to_id = {key: i+1 for i, (key, value) in enumerate(top_genera_to_stats.items())}
    write_out_extracted_sequences(top_genera_extracted_regions, top_genera_to_id, args.temp_dir)

    # simulate reads using ART from each extracted sequence file
    simulate_reads_from_extracted_sequences(top_genera_to_id, args.temp_dir)

    # calculate number of reads to extract from each genus
    genus_id_to_read_count = calculate_read_count_per_genus(top_genera_to_id, top_genera_to_stats, empty_genera, args.temp_dir)

    # gather the required # of reads from each genus and write to final file
    write_final_read_file(genus_id_to_read_count, args.temp_dir)

    # write the metadata that details for each ID what genus and how reads there are
    write_readset_metadata(top_genera_to_id, genus_id_to_read_count, args.temp_dir)

def process_reference_file(ref_path, genera_to_seqs):
    """ take reference sequences in SILVA file and parse into dictionary indexed by genera """

    def parse_taxonomy_lineage_to_genus(lineage):
        """
        e.g. Bacteria;Firmicutes;Bacilli;Bacillales;Planococcaceae;Planococcus;Planococcus citreus
             would return ...
             Bacteria;Firmicutes;Bacilli;Bacillales;Planococcaceae;Planococcus;
        """
        # lineage_split = lineage.split(";")
        # if "uncultured" not in lineage_split[-1] and "unidentified" not in lineage_split[-1]:
        #     assert lineage_split[-2] in lineage_split[-1], f"genus {lineage_split[-2]} not in species {lineage_split[-1]}"
        return ";".join(lineage.split(";")[:-1]) + ";"

    # go through each sequences in SILVA file
    with open(ref_path, "r") as ref_file:
        curr_seq = ""
        curr_lineage = ""
        included_count = 0; not_included_count = 0
        for line in ref_file:
            if line[0] == ">" and curr_seq != "":
                # handle previous seq
                if curr_lineage in genera_to_seqs:
                    genera_to_seqs[curr_lineage].append(curr_seq)
                    included_count += 1
                else:
                    not_included_count += 1

                # reset for next seq
                curr_seq = ""
                curr_lineage = parse_taxonomy_lineage_to_genus(" ".join(line.split()[1:]))
            elif line[0] == ">" and len(curr_seq) == 0:
                curr_lineage = parse_taxonomy_lineage_to_genus(" ".join(line.split()[1:]))
            else:
                curr_seq += line.strip()
        
        # process last seq
        if curr_lineage in genera_to_seqs:
            genera_to_seqs[curr_lineage].append(curr_seq)
            included_count += 1
        else:
            not_included_count += 1

    log_message("results", f"{included_count} sequences were specified" 
                           f" to genus level, {not_included_count} were not.\n")
    return genera_to_seqs

def process_primer_file(primer_file):
    """ process the primer file and return the forward and reverse primer sequences """

    # define tables to be used for complement and degenerate bases
    trans_tab = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C', '[': ']', ']': '['}
    deg_tab = {'M': '[AC]', 'Y': '[CU]', 'N': '[ACGU]', 'W': '[AU]', 'H': '[ACU]', 'V': '[ACG]', 'R': '[AG]'}

    # read file: first line is forward primer and second is reverse
    with open(primer_file, "r") as in_fd:
        for_primer = in_fd.readline().strip().upper().replace('T', 'U')
        rev_primer = in_fd.readline().strip().upper().replace('T', 'U')
    
    for_primer = "".join([deg_tab.get(base, base) for base in for_primer])
    rev_primer = "".join([deg_tab.get(base, base) for base in rev_primer])

    # reverse complement the reverse primer
    rev_primer = "".join([trans_tab.get(base, base) for base in rev_primer[::-1]])

    return for_primer, rev_primer

def process_biome_abundance_file(biome_abund_path):
    """ process the file generated by stats sub-command """

    # read in the file and store the top n most abundant genera
    with open(biome_abund_path, "r") as in_fd:
        all_lines = [line.strip().split("\t") for line in in_fd][1:]
        assert all([len(x) == 3 for x in all_lines]) # only 3 fields

        top_genera_to_stats = {line[0]: (int(line[1]), float(line[2])) for line in all_lines}
    return top_genera_to_stats

def extract_variable_regions(top_genera_extracted_regions, genera_to_seqs, for_primer, rev_primer):
    """ extract the variable regions using in-silico PCR for each top n genus """

    def mutate_sequence(seq, mutation_rate=0.02):
        dna_bases = ['A', 'C', 'G', 'T']
        return ''.join(base if random.random() > mutation_rate else random.choice(dna_bases) for base in seq)

    # set up the regex for looking primers (allow up to 3 subtitutions)
    for_pattern = regex.compile("(" + for_primer + "){s<=3}")
    rev_pattern = regex.compile("(" + rev_primer + "){s<=3}")

    # go through each genus and extract the variable region using primers
    total_length = 0; num_extractions = 0; empty_genera = []
    log_message("info", "extracting variable regions and randomly modifying 2% percent of those bases ...")

    for curr_genus in top_genera_extracted_regions.keys():
        assert curr_genus in genera_to_seqs, f"genus {curr_genus} not found in genera_to_seqs"
        
        num_extracted_for_current_genus = 0
        for curr_seq in genera_to_seqs[curr_genus]:
            # find forward and reverse primer matches
            for_matches = list(for_pattern.finditer(curr_seq))
            rev_matches = list(rev_pattern.finditer(curr_seq))

            # check if there is a single primer match for both
            if len(for_matches) == 1 and len(rev_matches) == 1:
                start = for_matches[0].end()
                end = rev_matches[0].start()

                # make sure end after start, and it seems to be a reasonable size
                if end > start and (end - start) > 200:
                    top_genera_extracted_regions[curr_genus].append(mutate_sequence(curr_seq[start:end].replace('U', 'T')))
                    total_length += (end - start)
                    num_extractions += 1
                    num_extracted_for_current_genus += 1
        
        # handle case where we not able to extract any variable regions
        if num_extracted_for_current_genus == 0:
            empty_genera.append(curr_genus)
    
    log_message("results", f"extracted {num_extractions} variable regions, total length: {total_length}, avg length: {total_length/num_extractions:.2f}\n")   
    return top_genera_extracted_regions, empty_genera

def write_out_extracted_sequences(top_genera_extracted_regions, top_genera_to_id, temp_dir):
    """ write out each genera's sequences as its own file to simulate reads """
    temp_dir = temp_dir if temp_dir[-1] == "/" else temp_dir + "/"

    log_message("info", "writing out extracted regions to file.\n")
    for genus, genus_num in top_genera_to_id.items():
        file_path = temp_dir + f"genus_{genus_num}_seqs.fna"

        if len(top_genera_extracted_regions[genus]) == 0:
            log_message("warning", f"no extracted regions for genus {genus}, skipping ...")
            continue

        with open(file_path, "w") as out_fd:
            for i, seq in enumerate(top_genera_extracted_regions[genus]):
                out_fd.write(f">genus_{genus_num}_seq_{i}\n{seq}\n")
    print()

def run_art_command(args):
    """ generates ART command-line and runs it """
    genus_num, temp_dir = args

    # only run if the file exists ... meaning we did extract sequences
    if not os.path.exists(f"{temp_dir}genus_{genus_num}_seqs.fna"):
        log_message("warning", f"no sequences found for genus {genus_num}, skipping ...")
        return

    command = f"art_illumina -ss MSv1 -amp -p -na -i {temp_dir}genus_{genus_num}_seqs.fna -l 250 -f 2000 -o {temp_dir}genus_{genus_num}_reads_"
    run_command(command)

def simulate_reads_from_extracted_sequences(top_genera_to_id, temp_dir):
    """ simulate reads from each genera in order to build final dataset """
    temp_dir = temp_dir if temp_dir[-1] == "/" else temp_dir + "/"
    log_message("info", "simulating reads from extracted sequences using ART ...\n")

    with Pool(20) as p:
        p.map(run_art_command, [(genus_num, temp_dir) for genus_num in top_genera_to_id.values()])
    print()

def calculate_read_count_per_genus(top_genera_to_id, top_genera_to_stats, empty_genera, temp_dir, total_read_count=200000):
    """ calculate how many reads we want from each genus """
 
    genus_id_to_read_count = {}; non_empty_genera = 0; empty_ids = []
    for genus, genus_id in top_genera_to_id.items():
        # add these extra checks I found scenarios where
        # no sequences with given primers are extracted
        # or ART simply does not simulate reads from a reference.
        if genus not in empty_genera and os.path.getsize(temp_dir + f"genus_{genus_id}_reads_1.fq") > 0:
            curr_ratio = top_genera_to_stats[genus][1]
            curr_read_count = int(curr_ratio * total_read_count)

            genus_id_to_read_count[genus_id] = curr_read_count
            total_read_count = total_read_count - curr_read_count
            non_empty_genera += 1
        else:
            genus_id_to_read_count[genus_id] = 0
            empty_ids.append(genus_id)
    
    # ... evenly distribute leftover reads
    batch_size = int(total_read_count/non_empty_genera)
    total_read_count -= batch_size * non_empty_genera

    genus_id_to_read_count = {genus_id: (read_count + batch_size if genus_id not in empty_ids else 0)  for genus_id, read_count in genus_id_to_read_count.items()}
    for genus_id in genus_id_to_read_count:
        if total_read_count > 0 and genus_id not in empty_ids:
            genus_id_to_read_count[genus_id] += 1
            total_read_count -= 1
    return genus_id_to_read_count

def write_final_read_file(genus_id_to_read_count, temp_dir):
    """ gather the required # of reads from each genus and write to final file """

    temp_dir = temp_dir if temp_dir[-1] == "/" else temp_dir + "/"
    log_message("info", "writing final read file ...")

    # write out the final read file
    with open(temp_dir + "final_reads_mate_1.fq", "w") as mate1, open(temp_dir + "final_reads_mate_2.fq", "w") as mate2:
        for genus_id, read_count in genus_id_to_read_count.items():
            total_reads_needed = genus_id_to_read_count[genus_id]
            curr_read_count = 1

            if total_reads_needed > 0:
                # gather all the reads for the current genus
                with open(temp_dir + f"genus_{genus_id}_reads_1.fq", "r") as in_mate1, open(temp_dir + f"genus_{genus_id}_reads_2.fq", "r") as in_mate2:
                    in_mate1_lines = [x.strip() for x in in_mate1.readlines()]
                    in_mate2_lines = [x.strip() for x in in_mate2.readlines()]

                    print(genus_id, total_reads_needed, len(in_mate1_lines), len(in_mate2_lines))

                    assert len(in_mate1_lines) == len(in_mate2_lines), "mate1 and mate2 files not same length"
                    assert len(in_mate1_lines) % 4 == 0, "mate1 is not a proper length"
                    assert len(in_mate1_lines) / 4 > total_reads_needed, "not enough reads in file"

                
                # randomly choose reads to include in final file
                selected_reads = random.sample(range(len(in_mate1_lines)//4), total_reads_needed)
                for i in selected_reads:
                    mate1.write(f"@genus_{genus_id}_read_{curr_read_count}_mate1\n{in_mate1_lines[i*4+1]}\n{in_mate1_lines[i*4+2]}\n{in_mate1_lines[i*4+3]}\n")
                    mate2.write(f"@genus_{genus_id}_read_{curr_read_count}_mate2\n{in_mate2_lines[i*4+1]}\n{in_mate2_lines[i*4+2]}\n{in_mate2_lines[i*4+3]}\n")
                    curr_read_count += 1

    log_message("info", f"final read file written to {temp_dir}final_reads_mate{{1,2}}.fna\n")

def write_readset_metadata(top_genera_to_id, genus_id_to_read_count, temp_dir):
    """ write the metadata file with taxonomy lineage, genus id, and read count """
    temp_dir = temp_dir if temp_dir[-1] == "/" else temp_dir + "/"
    log_message("info", "writing readset description file ...")

    with open(temp_dir + "seqtax.txt", "w") as out_fd:
        for genus, genus_id in top_genera_to_id.items():
            out_fd.write(f"{genus}\t{genus_id}\t{genus_id_to_read_count[genus_id]}\n")

    log_message("info", f"readset description file written to {temp_dir}seqtax.txt\n")