#!/usr/bin/env python3

import os
import sys
import subprocess
from datetime import datetime

# ANSI color codes for terminal output
IGREEN = '\033[0;92m'
IBLUE = '\033[0;94m'
IYELLOW = '\033[0;93m'
COLOR_OFF = '\033[0m'

# Paths to reference files (these will need to be customized for your environment)
ref_path_38 = '/path/to/GRCh38'
reference_38 = 'GRCh38.fa'
ref_path_19 = '/path/to/GRCh37'
reference_19 = 'GRCh37.fa'

project_path = '/Users/Patrick/Desktop/Python Projects/Projects'
rawfqs_path = '/Users/Patrick/Desktop/Python Projects/Projects/raw_fqs'  

   # Print message with a timestamp.
def echo_time(message: str):
    timestamp = f"{IGREEN}{datetime.now().strftime('%a %H:%M:%S')}{COLOR_OFF}"
    print(f"{timestamp} {message}")

    # Displays usage instructions for the script.
def display_usage():
    usage = """
    Usage: python3 script.py 
        -n | --name <project_name> 
        -s | --sequencing <data_type> 
        -r | --reference <reference_build>
    """
    print(usage)

    # Parses command-line arguments.
def check_arguments(args):
    if len(args) != 6:
        display_usage()
        sys.exit(1)
    
    arg_dict = {}
    i = 0
    while i < len(args):
        if args[i] in ('-n', '--name'):
            arg_dict['project_name'] = args[i + 1]
            i += 2
        elif args[i] in ('-s', '--sequencing'):
            data_type = args[i + 1].lower()
            if data_type in ['se', 'single', 'single-end']:
                arg_dict['data_type'] = 'SE'
            elif data_type in ['pe', 'paired', 'paired-end']:
                arg_dict['data_type'] = 'PE'
            else:
                echo_time("Sequencing type not recognized.")
                display_usage()
                sys.exit(1)
            i += 2
        elif args[i] in ('-r', '--reference'):
            reference_build = args[i + 1].lower()
            if reference_build in ['grch38', 'hg38']:
                arg_dict['ref_path'] = ref_path_38
                arg_dict['reference'] = reference_38
            elif reference_build in ['grch37', 'hg19']:
                arg_dict['ref_path'] = ref_path_19
                arg_dict['reference'] = reference_19
            else:
                echo_time("Reference build not recognized.")
                display_usage()
                sys.exit(1)
            i += 2
        else:
            echo_time(f"Unknown option: {args[i]}")
            display_usage()
            sys.exit(1)
    return arg_dict

    # Create project directories for output files.
def create_directories(project_name):
    project_dir = os.path.join(project_path, project_name)
    sub_dirs = ['fastp', 'bwa', 'bam', 'mapDamage', 'logs', 'statistics']
    if not os.path.exists(project_dir):
        echo_time(f"Creating directories for project {project_name}...")
        for sub_dir in sub_dirs:
            os.makedirs(os.path.join(project_dir, sub_dir), exist_ok=True)
    else:
        echo_time(f"Directories already exist for project {project_name}.")

    # Prompts user for the raw fastq directory and checks if it exists.
def check_fastq_directory():
    raw_fq = input(f"{IGREEN}Enter directory name of raw fq files: {COLOR_OFF}").strip()
    fq_path = os.path.join(rawfqs_path, raw_fq)
    if not os.path.isdir(fq_path):
        echo_time(f"{raw_fq} directory not found. Ensure the directory exists and try again.")
        sys.exit(1)
    else:
        echo_time("Creating list of fastqs to be processed...")
    return raw_fq

    # Creates a list of fastq files to be processed.
def create_fastq_list(project_name, raw_fq):
    fq_path = os.path.join(rawfqs_path, raw_fq)
    fastq_samples = sorted(
        set(["_".join(os.path.basename(f).split('_')[:2]) for f in os.listdir(fq_path) if f.endswith('.fastq.gz')])
    )

    fastq_list_path = os.path.join(project_path, project_name, f"{project_name}_fastqs.list")
    with open(fastq_list_path, 'w') as file:
        for sample in fastq_samples:
            file.write(f"{sample}\n")
    
    echo_time(f"Created fastq sample list: {IYELLOW}{fastq_list_path}{COLOR_OFF}")
    return fastq_samples

    # Creates a configuration file containing metadata for analysis. 
def create_config_file(project_name, arg_dict, raw_fq, fastq_samples):
    config_path = os.path.join(project_path, project_name, f"{project_name}.config")
    config_variables = {
        "project_name": project_name,
        "project_dir": os.path.join(project_path, project_name),
        "reference_build": arg_dict['reference'],
        "sequencing_type": arg_dict['data_type'],
        "rawfqs_path": os.path.join(rawfqs_path, raw_fq),
        "number_samples": len(fastq_samples),
        "reference_genome": os.path.join(arg_dict['ref_path'], arg_dict['reference']),
    }

    with open(config_path, 'w') as config_file:
        for key, value in config_variables.items():
            config_file.write(f'{key}="{value}"\n')
    
    echo_time(f"Created configuration file: {IYELLOW}{config_path}{COLOR_OFF}")

def main():
    echo_time("Starting aDNA analysis pipeline setup script")
    args = sys.argv[1:]
    arg_dict = check_arguments(args)

    project_name = arg_dict['project_name']
    data_type = arg_dict['data_type']
    reference_build = arg_dict['reference']
    create_directories(project_name)
    
    raw_fq = check_fastq_directory()
    fastq_samples = create_fastq_list(project_name, raw_fq)
    
    create_config_file(project_name, arg_dict, raw_fq, fastq_samples)
    echo_time("Setup script complete.")

if __name__ == '__main__':
    main()
