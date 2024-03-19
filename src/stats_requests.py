
# Name: stats_requests.py
# Description: This code contains methods related to stats sub-command
#              in the main script.

import biome_requests
import aiohttp
import asyncio
import pandas as pd
from taxonomy_path import TaxPath
from utils import log_message

def stats_main(args):
    # parse the SILVA taxonomy file and create a dictionary of genera paths
    fullpath_to_count, fingerprint_to_fullpath, full_lineage_dict = parse_taxonomy_file(args.taxonomy)

    # grab the file paths to studies for the given biome
    study_accessions = asyncio.run(biome_requests.grab_num_biome_studies(args.biome, print_results=False))
    log_message("info", f"found {len(study_accessions)} studies for this biome.\n")
    
    biome_file_paths = asyncio.run(biome_requests.grab_num_taxa_analyses(args.biome, study_accessions, print_results=False))
    log_message("info", f"found {len(biome_file_paths)} taxonomic analyses for this biome.\n")
    
    # now lets go through each analyses, and increment read counts for each genera
    fullpath_to_count = process_analyses_files(biome_file_paths, fullpath_to_count, fingerprint_to_fullpath)

    # write the output file with read counts
    write_output_file(fullpath_to_count, args.output)

def get_tax_lineage(line_split):
    """
    return: full taxonomic lineage from the SILVA taxonomy file

    The main special case we need to handle is when there are spaces
    in the taxonomic lineage, so we need to make sure we take care
    of that.

    e.g. Archaea;Crenarchaeota;Nitrososphaeria;Nitrosopumilales;Nitrosopumilaceae;Candidatus Nitrosopelagicus; 42941 genus 138
         should return ...
         Archaea;Crenarchaeota;Nitrososphaeria;Nitrosopumilales;Nitrosopumilaceae;Candidatus Nitrosopelagicus;
    """
    last_pos_of_traversal = max([i for i in range(len(line_split)) if ";" in line_split[i]])
    traversal = ' '.join(line_split[0:last_pos_of_traversal+1])
    return traversal

def get_tax_rank(line_split):
    """
    return: the taxonomic rank for a given line from SILVA file

    This is only slighly complicated by the fact that if you split
    every line by spaces, you could have a variable number. So you
    have go from the back BUT for some lines like Archaea; it looks like 
    this ...

    Archaea; 2 domain 

    Usually the is an additional field at the end that specifies the version.

    So the rule I will use is to scan for first non-number and return that 
    as the rank.
    """ 
    last_pos_of_traversal = max([i for i in range(len(line_split)) if ";" in line_split[i]])
    tax_rank = line_split[last_pos_of_traversal+2]
    return tax_rank 
    
def parse_taxonomy_file(taxonomy_file, ignore_eukaryotes=True):
    """ parse the SILVA taxonomy file """

    # read in the SILVA taxonomy file and build a dictionary mapping lineage to rank
    with open(taxonomy_file, 'r') as in_fd:
        all_lines = [line.strip().split() for line in in_fd]
        full_dict = {get_tax_lineage(x): get_tax_rank(x) for x in all_lines}
    assert len(all_lines) == len(full_dict), "issue with parsing the SILVA taxonomy file"

    # for each taxon specified to genus, build a TaxPath object
    genera_objs = [TaxPath(x, "silva", full_dict) for x in full_dict.keys() if full_dict[x] == "genus"]
    log_message("info", f"found {len(genera_objs)} genera in the SILVA taxonomy file.")

    ########################################################
    # important assumption: we will only focus on the
    # genera in Bacteria and Archaea domains.
    ########################################################
    if ignore_eukaryotes:
        log_message("info", f"important assumption: we will only focus genera in Bacteria and Archaea domains.")
        genera_objs = [x for x in genera_objs if x.domain_ in ["Bacteria", "Archaea"]]
        print()

    num_bacteria = len([x.domain_ for x in genera_objs if x.domain_ == "Bacteria"])
    num_archaea = len([x.domain_ for x in genera_objs if x.domain_ == "Archaea"])
    num_eukaryotes = len([x.domain_ for x in genera_objs if x.domain_ == "Eukaryota"])
    log_message("info", f"# genera in bacteria: {num_bacteria}, # genera in archaea: {num_archaea}, # genera in eukaryota: {num_eukaryotes}\n") 
    
    #########################################################
    # ultimately, we want to get a read count for each genera
    # so thats the first dictionary below. 
    #
    # however, we also know due to taxonomy discrepancies, we 
    # might not be able to figure out which genera reads belong 
    # to. So we can check if its "fingerprint" exists in the 
    # dictionary and use that to guide us a genus.
    #########################################################
    fullpath_to_count = {genera_objs[i].full_path: 0 for i in range(len(genera_objs))}
    
    # build a dictionary mapping fingerprint to full path
    fingerprint_to_fullpath_list = [(genera_objs[i].get_fingerprint(), genera_objs[i].full_path) for i in range(len(genera_objs))]
    fingerprint_to_fullpath_dict = {}; redundant_fingerprints = []

    # remove any redundant fingerprints
    for fingerprint, full_lineage in fingerprint_to_fullpath_list:
        if fingerprint not in fingerprint_to_fullpath_dict and fingerprint not in redundant_fingerprints:
            fingerprint_to_fullpath_dict[fingerprint] = full_lineage
        elif fingerprint in fingerprint_to_fullpath_dict:
            redundant_fingerprints.append(fingerprint)
            del fingerprint_to_fullpath_dict[fingerprint]
            log_message("warning", f"redundant fingerprint: {fingerprint}")
    assert len(fullpath_to_count) >= len(fingerprint_to_fullpath_dict), "issue with building the dictionaries"
    print()

    log_message("info", f"found {len(fullpath_to_count)} genera full lineages")
    log_message("info", f"found {len(fingerprint_to_fullpath_dict)} unique fingerprints\n")
    
    return fullpath_to_count, fingerprint_to_fullpath_dict, full_dict

def download_tsv(url):
    try:
        df = pd.read_csv(url, sep='\t')
        return df
    except Exception as e:
        print(f"\033[31mUnable to open file:\033[0m {url}. Error: {str(e)}")
        return None

def process_analyses_files(biome_file_paths, fullpath_to_count, fingerprint_to_fullpath):
    """ download each file, and increment the genera counts with read counts """
    direct_hits = 0; fp_hits = 0; no_hits = 0

    for curr_analysis_file in biome_file_paths:
        log_message("info", f"processing file: {curr_analysis_file}")

        # download file and sum up read counts across all samples
        df = download_tsv(curr_analysis_file)
        if df is None:
            continue
        df['total_read_count'] = df.iloc[:, 1:].sum(axis=1)

        # now lets go through each row and create a TaxPath object
        for index, row in df.iterrows():
            curr_classification = row.iloc[0]

            # check if classification is to genus level at least ...
            curr_obj = TaxPath(curr_classification, "ebi")
            if curr_obj.rank_level == "genus" and curr_obj.domain_ in ["Bacteria", "Archaea"]:
                curr_full_path = curr_obj.get_fullpath()
                curr_fp = curr_obj.get_fingerprint()

                if curr_full_path in fullpath_to_count:
                    fullpath_to_count[curr_full_path] += int(row['total_read_count'])
                    direct_hits += 1
                elif curr_fp in fingerprint_to_fullpath:
                    fullpath_to_count[fingerprint_to_fullpath[curr_fp]] += int(row['total_read_count'])
                    fp_hits += 1
                else:
                    no_hits += 1
                    #print("could not find a match for this classification: ", curr_full_path, curr_obj.rank_level)
    
    print()
    log_message("results", f"direct_hits: {direct_hits}, fingerprint_hits: {fp_hits}, no_hits: {no_hits}")
    return fullpath_to_count

def write_output_file(fullpath_to_count, output_path, top_n_genera=100):
    """ write the output file with read counts """

    # sort by read count, largest to smallest
    sorted_dict = dict(sorted(fullpath_to_count.items(), key=lambda item: item[1], reverse=True))
    total_read_count = sum(sorted_dict.values())

    # accumulate the top_n_genera results, or until the read count is zero
    output_results = []
    for i, (key, value) in enumerate(sorted_dict.items()):
        if value == 0 or len(output_results) == top_n_genera:
            break
        output_results.append((key, value, value/total_read_count))
    
    # throw error if we don't have at least 100 genera
    assert len(output_results) == top_n_genera, "issue with accumulating the top genera"

    # write the output file
    with open(output_path, 'w') as out_fd:
        out_fd.write("full_path,read_count,percent_total\n")
        for result in output_results:
            out_fd.write(f"{result[0]}\t{result[1]}\t{result[2]:.4f}\n")

