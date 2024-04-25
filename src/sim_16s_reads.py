#!/usr/bin/python3

# Name: sim_16s_reads.py
# Description: This script is meant to simulate 16S rRNA reads 
#              based on the most abundant genera in different 
#              types of biomes.

import argparse
from utils import valid_dir
from biome_requests import biome_main
from stats_requests import stats_main
from simulate import simulate_main

def parse_arguments():
    parser = argparse.ArgumentParser(description="simulate 16S rRNA reads, as well as general public data analysis.")
    sub_parser = parser.add_subparsers(dest='command')

    biome_parser = sub_parser.add_parser("biome", help="analyze biomes data in EBI database.")
    biome_parser.add_argument("--top-ten", dest="top_ten", action="store_true", help="report top ten most common biomes in EBI database.")
    biome_parser.add_argument("--grab-num-studies", dest="lineage", type=str, default="", help="find number of samples for certain biome")
    biome_parser.add_argument("--grab-num-taxa-analyses", dest="lineage_taxa", type=str, default="", help="find number of taxonomic analyses for certain biome")
    
    stats_parser = sub_parser.add_parser("stats", help="analyze biome samples in the context of SILVA taxonomy to identify most abundant genera.")
    stats_parser.add_argument("--biome", dest="biome", type=str, default="", help="biome lineage to analyze", required=True)  
    stats_parser.add_argument("--taxonomy", dest="taxonomy", type=str, default="", help="SILVA taxonomy file (*.txt)", required=True)  
    stats_parser.add_argument("--output", dest="output", type=str, default="output.tsv", help="output file for read counts for top_n genera")

    sim_parser = sub_parser.add_parser("simulate", help="simulate 16S rRNA reads based on most abundant genera in given biome.")
    sim_parser.add_argument("--biome-abundance", dest="biome_abund_path", type=str, default="", help="output from stats sub-command (*.tsv)", required=True)
    sim_parser.add_argument("--silva-ref", dest="silva_ref_path", type=str, default="", help="SILVA reference file (*.fasta)", required=True)
    sim_parser.add_argument("--silva-taxonomy", dest="silva_tax_path", type=str, default="", help="SILVA taxonomy file (*.txt)", required=True)
    sim_parser.add_argument("--num-reads", dest="num_reads", type=int, default=200000, help="number of reads to simulate", required=True)
    sim_parser.add_argument("--output-name", dest="output_name", type=str, help="output file name prefix for simulated reads", required=True)
    sim_parser.add_argument("--primers", dest="primer_file", type=str, default="", help="file containing forward and reverse primer sequences", required=True)
    sim_parser.add_argument("--temp-dir", dest="temp_dir", type=valid_dir, default="", help="temporary directory for intermediate files", required=True)

    args = parser.parse_args()

    # handle situations where help should be printed
    if args.command is None:
        parser.print_help()
        exit(1)
    elif args.command == "biome" and not (args.top_ten or args.lineage or args.lineage_taxa):
        biome_parser.print_help()
        exit(1)
    return args

if __name__ == "__main__":
    args = parse_arguments(); print()

    if args.command == "biome":
        biome_main(args)
    elif args.command == "stats":
        stats_main(args)
    elif args.command == "simulate":
        simulate_main(args)
    print()