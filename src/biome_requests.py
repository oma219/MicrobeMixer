
# Name: biome_requests.py
# Description: This code contains methods related to biome sub-command
#              in the main script.

import time
import multiprocessing as mp
import aiohttp
import asyncio

from utils import log_message, get_request, get_request_async
import multiprocessing as mp

def biome_main(args):
    if args.top_ten:
        get_top_ten_biomes()
    elif args.lineage != "":
        asyncio.run(grab_num_biome_studies(args.lineage))
    elif args.lineage_taxa != "":
        study_accessions = asyncio.run(grab_num_biome_studies(args.lineage_taxa, print_results=False))
        log_message("info", f"found {len(study_accessions)} studies for this biome.\n")
        asyncio.run(grab_num_taxa_analyses(args.lineage_taxa, study_accessions))

def get_top_ten_biomes():
    """ report the top ten most common biomes in EBI database. """

    request_url = "https://www.ebi.ac.uk/metagenomics/api/v1/biomes/top10?format=json"
    log_message("info", f"requesting data from {request_url}\n")
    
    # get results from EBI API
    json_response = get_request(request_url)

    # iterate through API response
    for i, biome_json in enumerate(json_response['data']):
        biome_name = biome_json['attributes']['biome-name']
        num_samples = biome_json['attributes']['samples-count']
        lineage = biome_json['attributes']['lineage']
        
        log_message("results", f"biome #{i}, name: {biome_name:17} (samples={num_samples}), lineage: {lineage:40s}")

async def grab_num_biome_studies(lineage, print_results=True):
    """ identify the number of samples for certain biome in EBI database. """
    
    log_message("info", f"requesting all studies from EBI API for biome: {lineage}\n")

    # create the request url for this task
    lineage = lineage.replace(':', '%3A'); page = 1
    base_url = f"https://www.ebi.ac.uk/metagenomics/api/v1/biomes/{lineage}%20/studies?format=json&page="

    async with aiohttp.ClientSession() as session:
        # get the first page to find out how many pages there are
        json_response = await get_request_async(session, base_url + str(page))
        total_pages = json_response['meta']['pagination']['pages']

        # start all the requests
        tasks = [get_request_async(session, base_url + str(i)) for i in range(2, total_pages + 1)]
        pages = await asyncio.gather(*tasks)

    # process the results
    study_list = json_response['data'] + [item for page in pages for item in page['data']]

    # get number of samples per study 
    accession_names = [study['attributes']['accession'] for study in study_list if study['attributes']['accession'] is not None]
    sample_count = [int(study['attributes']['samples-count']) for study in study_list if study['attributes']['samples-count'] is not None]
    
    if print_results:
        log_message("results", f"total number of samples: {sum(sample_count)}, number of studies: {len(sample_count)}")
        log_message("results", f"first ten accession names: {accession_names[:10]}\n")

    return accession_names

async def grab_num_taxa_analyses(lineage_taxa, study_accessions, print_results=True):
    """ identify how many studies have SSU taxa analyses for the latest version (v5.0) """

    log_message("info", f"iterating through all studies and checking for analyses file/version # ...\n")
    selected_files = []

    async with aiohttp.ClientSession() as session:
        # Process the studies in batches to avoid overloading the server
        batch_size = 10
        for i in range(0, len(study_accessions), batch_size):
            batch = study_accessions[i:i+batch_size]
            tasks = [process_study(session, curr_study) for curr_study in batch]
            try:
                results = await asyncio.gather(*tasks)
                selected_files += results
            except Exception as e:
                print(f"An error occurred: {e}")
            log_message("info", f"processed {i+batch_size} studies...")
    print()
    
    selected_files = [file for file in selected_files if file is not None]

    if print_results:
        log_message("results", f"total number of studies with SSU taxa analyses: {len(selected_files)}")
        log_message("results", f"first ten download links: {selected_files[:10]}\n")
    
    return selected_files

async def process_study(session, curr_study):
    """ looks through a specific study's downloadable files for SSU Analysis """
    base_url = f"https://www.ebi.ac.uk/metagenomics/api/v1/studies/{curr_study}/downloads?page="
    page = 1

    # get the first page to find out how many pages there are
    json_response = await get_request_async(session, base_url + str(page))
    total_pages = json_response['meta']['pagination']['pages']

    # start all the requests
    tasks = [get_request_async(session, base_url + str(i)) for i in range(2, total_pages + 1)]
    pages = await asyncio.gather(*tasks)
            
    # process the results
    data_list = json_response['data'] + [item for page in pages for item in page['data']]

    # check if SSU taxa analysis is present and it is the latest version
    for download_file in data_list:
        if download_file['attributes']['description']['label'] == "Taxonomic assignments SSU":
            if download_file['relationships']['pipeline']['data']['id'] == "5.0":
                return download_file['links']['self']

