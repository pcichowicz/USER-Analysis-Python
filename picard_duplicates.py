import os
import subprocess
from datetime import datetime

def echo_time(message):
    """Prints a message with a timestamp, similar to echo_time function in the bash script."""
    timestamp = datetime.now().strftime('%a %H:%M:%S')
    base = os.path.basename(__file__).replace(".py", "")
    print(f"\033[0;32m{timestamp}\033[0m \033[0;34m{base}\033[0m {message}")

# 1. Load configuration variables from the config file
project_name = "your_project_name_here"  # Replace with project name argument, like sys.argv[1]
project_dir = f"/Users/Patrick/aDNA/Project/{project_name}"
config_path = os.path.join(project_dir, f"{project_name}.config")

# Parse the config file and store key-value pairs in a dictionary
config = {}
with open(config_path) as config_file:
    for line in config_file:
        line = line.strip()
        if not line or line.startswith("#") or '=' not in line:
            continue  # Skip comments and invalid lines
        key, value = line.split("=", 1)
        config[key] = value.strip('"')

# Source required variables from config
logs_path = config.get('logs_path', f"{project_dir}/logs")
bam_path = config.get('bam_path', f"{project_dir}/bam")
number_samples = int(config.get('number_samples', 0))

# 2. Start the process for marking duplicates
echo_time("Running Picard MarkDuplicates")

markduplicates_counter = 1
fastq_list_path = os.path.join(project_dir, f"{project_name}_fastqs.list")

with open(fastq_list_path, 'r') as fastq_list_file:
    for mapped in fastq_list_file:
        mapped = mapped.strip()
        name = os.path.basename(mapped).replace(".fastq.gz", "")

        I = os.path.join(bam_path, f"{name}_mapped.bam")
        O = os.path.join(bam_path, f"{name}_mapped-md.bam")
        opt_duplicates = "-OPTICAL_DUPLICATE_PIXEL_DISTANCE 12000"
        remove_duplicates = "-REMOVE_DUPLICATES false"  # Change to 'true' if needed
        metrics_files = f"-METRICS_FILE {os.path.join(bam_path, f'{name}_mapped-md.bam.metrics')}"
        tag_policy = "-TAGGING_POLICY All"
        validation = "-VALIDATION_STRINGENCY LENIENT"

        echo_time(f"> Sample {markduplicates_counter} out of {number_samples} - {name}")

        # Picard command
        picard_command = [
            "java", "-jar", 
            "/usr/local/Cellar/picard/picard.jar", 
            "MarkDuplicates", 
            "-I", I, 
            "-O", O, 
            opt_duplicates, 
            remove_duplicates, 
            metrics_files, 
            tag_policy, 
            validation
        ]

        command_str = ' '.join(picard_command)
        echo_time(f"Running command: {command_str}")
        
        # Run the Picard MarkDuplicates command
        try:
            subprocess.run(picard_command, check=True)
        except subprocess.CalledProcessError as e:
            echo_time(f"Error running Picard MarkDuplicates for sample {name}: {e}")
            continue  # Skip to next sample

        # 3. Index the marked duplicate BAM file
        samtools_command = ["samtools", "index", "-b", O]
        echo_time(f"Running command: {' '.join(samtools_command)}")
        
        try:
            subprocess.run(samtools_command, check=True)
        except subprocess.CalledProcessError as e:
            echo_time(f"Error running samtools index for {O}: {e}")
        
        markduplicates_counter += 1

echo_time("Finished with marking and removing duplicate reads.")
