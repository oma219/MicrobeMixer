## 🦠 MicrobeMixer: 16S rRNA Read Dataset Simulator

## Description:

This code provides utilities to **simulate realistic 16S rRNA read datasets from different environments** based on public data. The code will generate datasets using
`EBI Mgnify` service to query public data, and summarize the most common bacteria present in the specified biome. Followed by extracting variables with
user-given primers and simulating reads from those variable regions.

## Getting Started:
In order to start running the code, you just have to clone the repository and install a couple of packages. These packages include `regex`, `aiohttp` and `asyncio` 
and they can be installed using either `pip` or `conda/mamba`.

For running the sample datasets, you will need to download the `SILVA` reference file which also shown below in the example commands.

```base
git clone https://github.com/oma219/MicrobeMixer.git

cd MicrobeMixer/data
wget https://ftp.arb-silva.de/release_138.1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
gzip -d SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
```

## How to Simulate a Read Set?

To simulate reads, there are 2 main steps. First, specify a biome and summarize the most common bacteria. Second, simulate reads from the variable
regions in the rRNA for the most common bacteria.

#### Step 1: Find the most common bacteria in a certain biome

For the example, let's say we want to simulate a realistic `Human Gut` 16S rRNA read dataset. We would start with finding the most common bacteria 
in human gut samples in EBI by using the command below.

```
python3 sim_16s_reads.py stats 
                         --biome 'root:Host-associated:Human:Digestive system'  
                         --taxonomy ../data/tax_slv_ssu_138.1.txt
```
The string `root:Host-associated:Human:Digestive system` represents the EBI lineage for Human Gut Biome. In order to determine you can look at the `Other Utilities`
section below. The output of this command will be a file in the current directory called `output.tsv` which will be used in Step 2.

#### Step 2: Simulate a paired-end readset for a certain biome

Now lets simulate reads for this particular biome using this command:
```
python3 sim_16s_reads.py simulate
                         --biome-abundance output.tsv 
                         --silva-ref ../data/SILVA_138.1_SSURef_NR99_tax_silva.fasta 
                         --silva-taxonomy ../data/tax_slv_ssu_138.1.txt 
                         --primers ../data/V3_V4_primers.txt  
                         --temp-dir <OUTPUT_DIR> # Fill this in!
```
The output of this command will be two FASTQ files (`final_reads_mate_{1,2}.fq`) and metadata file (`seqtax.txt`) in the specified output directory. Within the
`data/`, there are other primer files that can be used or a user-specified example that follows the same format.

## Other Utilities:

- Find the top 10 most common biome in EBI (use this to find common biome lineages):
  ```
  python3 sim_16s_reads.py biome --top-ten
  ```

- Find the number of taxonomic analyses for studies in a certain biome in EBI:
  ```
  python3 sim_16s_reads.py biome --grab-num-taxa-analyses 'root:Host-associated:Human:Digestive system’
  ```
