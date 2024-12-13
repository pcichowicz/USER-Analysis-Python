import os
import subprocess
import sys
import datetime

def echo_time(message):
    """Prints a timestamped message to the console."""
    timestamp = datetime.datetime.now().strftime("%a %H:%M:%S")
    print(f"{timestamp} [INFO] {message}")

def run_command(command, log_file=None):
    """Run a shell command, log output to the console and optionally to a log file."""
    try:
        echo_time(f"Running command: {command}")
        result = subprocess.run(command, shell=True, check=True, text=True, capture_output=True)
        output = result.stdout.strip()
        
        if log_file:
            with open(log_file, 'a') as log:
                log.write(f"{output}\n")
        
        return output
    except subprocess.CalledProcessError as e:
        error_message = f"Command failed: {command}\nError: {e.stderr.strip()}"
        if log_file:
            with open(log_file, 'a') as log:
                log.write(f"{error_message}\n")
        print(f"[ERROR] {error_message}")
        sys.exit(1)

def load_config(config_path):
    """Load the config file into a dictionary (like 'source' in Bash)."""
    config = {}
    with open(config_path, 'r') as file:
        for line in file:
            if '=' in line:
                key, value = line.strip().split('=', 1)
                config[key.strip()] = value.strip()
    return config

def main(project_name):
    project_dir = f"/Users/Patrick/aDNA/Project/{project_name}"
    config_path = os.path.join(project_dir, f"{project_name}.config")
    log_path = os.path.join(project_dir, "logs")
    fastq_list_path = os.path.join(project_dir, f"{project_name}_fastqs.list")

    # Load configuration from the config file
    config = load_config(config_path)
    bam_path = config.get('bam_path', '')
    mapDamage_path = config.get('mapDamage_path', '')
    reference_genome = config.get('reference_genome', '')

    # Ensure log directory exists
    os.makedirs(log_path, exist_ok=True)
    log_file = os.path.join(log_path, f"{project_name}_mapDamage.log")

    echo_time("Running mapDamage script ...")

    try:
        with open(fastq_list_path, 'r') as sample_list:
            mapDam_counter = 1
            samples = sample_list.readlines()
            number_samples = len(samples)

            for sample in samples:
                name = os.path.basename(sample.strip())
                bam_input = os.path.join(bam_path, f"{name}_mapped-md.bam")
                folder = os.path.join(mapDamage_path, f"{name}_mapdamage")

                echo_time(f"> Sample {mapDam_counter} out of {number_samples}; Processing {name}")
                log_message = f"> Sample {mapDam_counter} out of {number_samples}; Processing {name}\n"
                
                # Log to the log file as well
                with open(log_file, 'a') as log:
                    log.write(log_message)

                map_damage_command = f"mapDamage -Q 30 -d {folder} -i {bam_input} -r {reference_genome}"
                echo_time(f"@ {map_damage_command}")
                
                # Log the command as well
                with open(log_file, 'a') as log:
                    log.write(f"@ {map_damage_command}\n")
                
                # Run the mapDamage command
                run_command(map_damage_command, log_file)
                
                mapDam_counter += 1
    except FileNotFoundError as e:
        print(f"[ERROR] File not found: {e}")
    except Exception as e:
        print(f"[ERROR] An error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 script.py <project_name>")
        sys.exit(1)
    project_name = sys.argv[1]
    main(project_name)
