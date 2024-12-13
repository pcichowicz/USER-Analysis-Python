#!/usr/bin/env python3

import sys
import os
import subprocess
from datetime import datetime

# ANSI color codes for terminal output
IGREEN = '\033[0;92m'
IBLUE = '\033[0;94m'
COLOR_OFF = '\033[0m'

# Script names mapped to their file paths
SCRIPTS = {
    'fastp': './fastp.sh',
    'bwa_aln': './bwa_aln.sh',
    'mark_duplicates': './mark_duplicates.sh',
    'mapdamage': './mapdamage.sh'
}

def echo_time(message: str):
    """
    Prints a message with a timestamp and the script name.
    """
    base = f"{IBLUE}{os.path.basename(__file__)}{COLOR_OFF}"
    timestamp = f"{IGREEN}{datetime.now().strftime('%a %H:%M:%S')}{COLOR_OFF}"
    print(f"{timestamp} {base} {message}")

def get_project_name():
    """
    Prompts the user for a project name and ensures it's not empty.
    """
    project_name = input(echo_time('Enter project name for config file: '))
    if not project_name.strip():
        echo_time('Project name cannot be empty. Exiting.')
        sys.exit(1)
    return project_name.strip()

def validate_scripts(arguments):
    """
    Validate the scripts passed as command-line arguments.
    Returns a list of valid script names and a list of invalid ones.
    """
    valid_scripts = []
    invalid_scripts = []

    for arg in arguments:
        if arg in SCRIPTS:
            valid_scripts.append(arg)
        else:
            invalid_scripts.append(arg)

    if invalid_scripts:
        echo_time(f"Invalid script names provided: {', '.join(invalid_scripts)}")
        echo_time(f"Valid options: {', '.join(SCRIPTS.keys())}")
    
    return valid_scripts

def run_scripts(valid_scripts, config_file):
    """
    Run the specified scripts, passing the config file as an argument.
    """
    for script in valid_scripts:
        script_path = SCRIPTS[script]
        echo_time(f"Running {script_path} script...")
        try:
            result = subprocess.run([script_path, config_file], check=True, text=True, capture_output=True)
            echo_time(f"Output from {script_path}:\n{result.stdout}")
        except subprocess.CalledProcessError as e:
            echo_time(f"Error while running {script_path}: {e.stderr}")
        except FileNotFoundError:
            echo_time(f"Script not found: {script_path}. Make sure it exists and is executable.")

def main():
    """
    Main function to handle script execution.
    """
    if len(sys.argv) < 2:
        echo_time("Please provide at least one script name as a parameter.")
        echo_time("Valid options: fastp | bwa_aln | mark_duplicates | mapdamage")
        sys.exit(1)

    # Get the project name from user input
    project_name = get_project_name()

    # Validate and get the list of scripts to run
    arguments = sys.argv[1:]
    valid_scripts = validate_scripts(arguments)
    
    if not valid_scripts:
        echo_time("No valid scripts provided. Exiting.")
        sys.exit(1)
    
    # Run the selected scripts
    run_scripts(valid_scripts, project_name)

    echo_time("Running of all script(s) complete. Log files are located in the project directory.")

if __name__ == '__main__':
    main()
