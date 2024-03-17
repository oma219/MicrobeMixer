# Name: utils.py
# Description: This code contains general methods
#              to be used for the entire project.

import requests
import aiohttp
import asyncio
import os
import argparse
import subprocess
from requests.exceptions import HTTPError

def log_message(level="info", message=""):
    if level == "error":
        print(f"\033[31m[log::{level}]\033[0m {message}")
    elif level == "warning":
        print(f"\033[33m[log::{level}]\033[0m {message}")
    else:
        print(f"\033[32m[log::{level}]\033[0m {message}")

def valid_dir(path):
    if os.path.isdir(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} is not a valid directory")

def get_request(request_url):
    """ performs a get request and returns results """    
    try:
        response = requests.get(request_url)
        response.raise_for_status()
        json_response = response.json()
    except HTTPError as http_err:
        log_message("error", f'HTTP error occurred: {http_err}'); exit(1)
    except Exception as err:
        log_message("error", f'Other error occurred: {err}'); exit(1)
    return json_response

async def get_request_async(session, url, retries=3):
    for i in range(retries + 1):
        try:
            async with session.get(url) as resp:
                resp.raise_for_status()  # Raises an exception for non-200 status codes
                return await resp.json(content_type=None)
        except aiohttp.ClientResponseError as e:
            if e.status == 502 and i < retries:  # Bad Gateway error
                print(f"\033[31m[Bad Gateway Error]\033[0m url: {url}, retrying...\n")
                await asyncio.sleep(1)  # Wait for 1 second before retrying
            else:
                print(f"Request to {url} returned an error: {e}")
                exit(1) # Exit if we've exhausted our retries or if it's a different error
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            exit(1)

def run_command(command):
    """ run shell command and check for errors """
    result = subprocess.run(command, shell=True, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        print(f"Command '{command}' failed with error code {result.returncode}")
        print(f"Output: {result.stdout.decode()}")
        print(f"Error: {result.stderr.decode()}")
    else:
        log_message("shell", f"Command '{command}' succeeded")
