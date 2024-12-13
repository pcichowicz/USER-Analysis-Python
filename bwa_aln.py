import os
import subprocess
from datetime import datetime

# Function to print timestamps with echo-like formatting
def echo_time(message):
    timestamp = datetime.now().strftime("%a %H:%M:%S")
    base = os.path.basename(__file__).replace(".py", "")
    print(f"\033[1;34m{timestamp}\033[0m \033[1;34m{base}\033[0m {message}")

# Paths and configurations
project_name = "your_project_name"  # Replace with your project name or dynamically fetch from input args
project_dir = f"/Users/Patrick/aDNA/Project/{project_name}"


# Load configuration from a file
config_path = os.path.join(project_dir, f"{project_name}.config")
config = {}

# with open(config_path) as config_file:
#     for line in config_file:
#         key, value = line.strip().split("=")
#         config[key] = value.strip('"')

with open(config_path) as config_file:
    for line in config_file:
        line = line.strip()
        # Skip empty lines and comments
        if not line or line.startswith("#"):
            continue
        # Make sure line contains an '='
        if '=' not in line:
            continue
        try:
            key, value = line.split("=", 1)  # Split only at the first '='
            config[key] = value.strip('"')
        except ValueError:
            print(f"Warning: Malformed line in config file: {line}")
# # Load color variables (if needed)
# color_path = "/Users/Patrick/aDNA/temp/colors.txt"
# colors = {}
# with open(color_path) as color_file:
#     for line in color_file:
#         if "=" in line:
#             key, value = line.strip().split("=")
#             colors[key] = value

# Assign configuration variables
fastp_path = config.get("fastp_path")
bwa_path = config.get("bwa_path")
bam_path = config.get("bam_path")
reference_genome = config.get("reference_genome")
logs_path = config.get("logs_path")
number_samples = int(config.get("number_samples", 0))
sequencing_platform = config.get("sequencing_platform", "").upper()

# Log file paths
bwa_aln_log = os.path.join(logs_path, f"{project_name}_bwa-aln.log")
bam_log = os.path.join(logs_path, f"{project_name}_bam.log")

echo_time(f"Running alignment script - BWA aln & samtools\n")
aln_counter = 1

# Read fastq list from project directory
fastq_list_path = os.path.join(project_dir, f"{project_name}_fastqs.list")
with open(fastq_list_path) as fastq_file:
    fastq_files = fastq_file.readlines()

for trimmed in fastq_files:
    trimmed = trimmed.strip()
    name = os.path.basename(trimmed).replace(".fastq.gz", "")

    merged = f"{fastp_path}/{name}_trimmed_merged.fastq.gz"
    R1 = f"{fastp_path}/{name}_trimmed_out_1.fastq.gz"
    R2 = f"{fastp_path}/{name}_trimmed_out_2.fastq.gz"
    S1 = f"{bwa_path}/{name}_trimmed_out_1.fastq.gz.sai"
    S2 = f"{bwa_path}/{name}_trimmed_out_2.fastq.gz.sai"
    M1 = f"{bwa_path}/{name}_trimmed_merged.fastq.gz.sai"

    # Sample metadata for read group (RG) tags
    ID = f"{name}_merged"
    ID2 = f"{name}_paired"
    SM = name[:7]
    CN = "CGG"
    PL = sequencing_platform.upper()
    LB = name.split('_')[3] if len(name.split('_')) > 3 else 'unknown'
    DS = "aDNA USER Reduction project"

    # Merged reads alignment
    echo_time(f"> Sample {aln_counter} out of {number_samples}; merged reads - {name}")
    bwa_aln_command = f"bwa aln {reference_genome} {merged} > {M1}"
    echo_time(bwa_aln_command)
    subprocess.run(bwa_aln_command, shell=True)

    # Paired-end reads alignment
    echo_time(f"> Sample {aln_counter} out of {number_samples}; paired-end reads 1 & 2 - {name}")
    bwa_aln_command_r1 = f"bwa aln {reference_genome} {R1} > {S1}"
    bwa_aln_command_r2 = f"bwa aln {reference_genome} {R2} > {S2}"
    echo_time(bwa_aln_command_r1)
    subprocess.run(bwa_aln_command_r1, shell=True)
    echo_time(bwa_aln_command_r2)
    subprocess.run(bwa_aln_command_r2, shell=True)

    # BWA samse command
    echo_time(f"> Sample {aln_counter} out of {number_samples}; samse - {name}")
    bwa_samse_command = (
        f"bwa samse -r '@RG\\tID:{ID}\\tSM:{SM}\\tCN:{CN}\\tPL:{PL}\\tLB:{LB}\\tDS:{DS}' "
        f"{reference_genome} {M1} {merged} | samtools sort -o {bam_path}/{name}_merged.bam"
    )
    echo_time(bwa_samse_command)
    subprocess.run(bwa_samse_command, shell=True)

    # BWA sampe command
    echo_time(f"> Sample {aln_counter} out of {number_samples}; sampe - {name}")
    bwa_sampe_command = (
        f"bwa sampe -r '@RG\\tID:{ID2}\\tSM:{SM}\\tCN:{CN}\\tPL:{PL}\\tLB:{LB}\\tDS:{DS}' "
        f"{reference_genome} {S1} {S2} {R1} {R2} | samtools sort -o {bam_path}/{name}_paired-reads.bam"
    )
    echo_time(bwa_sampe_command)
    subprocess.run(bwa_sampe_command, shell=True)

    # Merge BAM files
    echo_time(f"> Sample {aln_counter} out of {number_samples}; samtools merge - {name}")
    merge_command = f"samtools merge {bam_path}/{name}.bam {bam_path}/{name}*.bam"
    echo_time(merge_command)
    subprocess.run(merge_command, shell=True)

    # Extract mapped reads
    echo_time(f"> Sample {aln_counter} out of {number_samples}; extract mapped reads - {name}")
    extract_command = f"samtools view -b -F 4 -o {bam_path}/{name}_mapped.bam {bam_path}/{name}.bam"
    echo_time(extract_command)
    subprocess.run(extract_command, shell=True)

    # Index BAM file
    echo_time(f"> Sample {aln_counter} out of {number_samples}; index BAM file - {name}")
    index_command = f"samtools index -b {bam_path}/{name}_mapped.bam"
    echo_time(index_command)
    subprocess.run(index_command, shell=True)

    aln_counter += 1

echo_time("Finished read alignments and creating BAM files")
