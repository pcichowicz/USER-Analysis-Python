import os
import subprocess
from datetime import datetime

def echo_time(message):
    """Prints a message with a timestamp, similar to echo_time function in the bash script."""
    timestamp = datetime.now().strftime('%a %H:%M:%S')
    base = os.path.basename(__file__).replace(".py", "")
    print(f"\033[0;32m{timestamp}\033[0m \033[0;34m{base}\033[0m {message}")

# 1. Load configuration variables from the config file
project_name = "your_project_name_here"  # Replace this with the project name as needed
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

# Source required config values
logs_path = config.get('logs_path', f"{project_dir}/logs")
fastp_path = config.get('fastp_path', f"{project_dir}/fastp")
sequencing_type = config.get('sequencing_type', 'PE')
number_samples = int(config.get('number_samples', 0))
fastp = config.get('fastp', 'fastp')  # Assuming fastp is in PATH, otherwise add full path

# 2. Start of the processing logic
echo_time(f"Running adapter trimming and quality control script for project {project_name}")

fastq_list_path = os.path.join(project_dir, f"{project_name}_fastqs.list")
fastq_counter = 1

with open(fastq_list_path, 'r') as fastq_list_file:
    for samples in fastq_list_file:
        samples = samples.strip()
        name = os.path.basename(samples).replace(".fastq.gz", "")

        echo_time(f"> Sample {fastq_counter} out of {number_samples} - {name}")

        # 3. Define fastp command options
        in1 = f"--in1 {samples}_S1_L004_R1_001.fastq.gz"
        in2 = f"--in2 {samples}_S1_L004_R2_001.fastq.gz"
        out1 = f"--out1 {os.path.join(fastp_path, f'{name}_trimmed_out_1.fastq.gz')}"
        out2 = f"--out2 {os.path.join(fastp_path, f'{name}_trimmed_out_2.fastq.gz')}"
        unpair_1 = f"--unpaired1 {os.path.join(fastp_path, f'{name}_trimmed_unpaired.fastq.gz')}"
        unpair_2 = f"--unpaired2 {os.path.join(fastp_path, f'{name}_trimmed_unpaired.fastq.gz')}"
        merge = "--merge"
        merge_out = f"--merged_out {os.path.join(fastp_path, f'{name}_trimmed_merged.fastq.gz')}"
        failed_out = f"--failed_out {os.path.join(fastp_path, f'{name}_trimmed_failed_out.fastq.gz')}"
        pe_adapter = "--detect_adapter_for_pe"
        dont_eval = "--dont_eval_duplication"
        length_required = "--length_required 30"
        quality_phred = "--qualified_quality_phred 30"
        report_title = f"--report_title '{name}'"
        html = f"--html {os.path.join(fastp_path, f'{name}.html')}"
        json = f"--json {os.path.join(fastp_path, f'{name}.json')}"
        threads = "--thread 4"  # Assuming 4 threads as an example

        # 4. Execute fastp command for SE or PE sequencing
        if sequencing_type == "SE":
            command = [
                fastp, 
                length_required, quality_phred, in1, out1, 
                failed_out, dont_eval, report_title, html, json
            ]
        else:  # Paired-end
            command = [
                fastp, 
                length_required, quality_phred, in1, in2, out1, out2, 
                unpair_1, unpair_2, merge, merge_out, 
                failed_out, pe_adapter, dont_eval, 
                report_title, html, json
            ]

        command_str = ' '.join(command)
        echo_time(f"@ {command_str}\n")
        
        # Execute the fastp command
        try:
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as e:
            echo_time(f"Error running fastp for sample {name}: {e}")
        
        fastq_counter += 1

echo_time("fastp trimming and quality control finished.")
