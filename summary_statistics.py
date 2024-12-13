import os
import subprocess
import csv
import json
import sys

def run_command(command):
    """Run a shell command and return the output as a string."""
    try:
        result = subprocess.run(command, shell=True, check=True, text=True, capture_output=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {command}\n{e.stderr}")
        sys.exit(1)

def calculate_proportion(numerator, denominator):
    """Calculate the proportion, handle division by zero."""
    try:
        return round(numerator / denominator, 5) if denominator > 0 else 0
    except ZeroDivisionError:
        return 0

def main(project_name):
    project_dir = f"/Users/Patrick/aDNA/Project/{project_name}"
    config_path = f"{project_dir}/{project_name}.config"

    # Load configuration file (source equivalent)
    with open(config_path) as f:
        config = dict(line.strip().split('=') for line in f if '=' in line)
    
    # Paths from the config
    bam_path = config.get('bam_path', '')

    # Write the CSV header
    header = [
        "Sample ID", "Sample Age", "Library", "Treatment", 
        "Total Reads", "Trimmed Reads", "Mapped Reads", 
        "Cluster Duplicates", "PCR Duplicates", "Unique Mapped", 
        "pTrimmed Reads", "pMapped Reads", "Clonality", 
        "pCluster Duplicates", "pPCR duplicates", "Endogenous", 
        "pUnique Mapped", "Efficiency", "Read Length", 
        "xDepth", "yDepth", "mDepth", "autDepth"
    ]

    output_path = os.path.join(project_dir, 'stat.csv')

    with open(output_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(header)

        with open(f"{project_dir}/{project_name}_fastqs.list", 'r') as sample_list:
            for sample_path in sample_list:
                sample_path = sample_path.strip()
                name = os.path.basename(sample_path)

                sample_id = name.split('_')[0]
                sample_age = name.split('_')[1]
                library = name.split('_')[3]
                treatment = name.split('_')[5]

                # Total genome length
                bam_file = os.path.join(bam_path, f"{name}_mapped-md.bam")
                genome_length_cmd = f"samtools view -H {bam_file} | awk -vFS=: '/^@SQ/ {{sum+=$3}} END {{print sum}}'"
                genome_length = int(run_command(genome_length_cmd))

                # Calculate depth using samtools depth and an AWK script for X, Y, M, and autosomes
                depth_cmd = f"samtools depth -q30 -Q20 -a {bam_file} | head -1000 | awk -f sexDetermination.awk"
                depth = list(map(float, run_command(depth_cmd).split()))

                # Calculate total, trimmed, and mapped reads
                total_reads = int(run_command("jq '.summary.before_filtering.total_reads' fastp.json"))
                trimmed_reads = int(run_command(f"samtools view -c {bam_path}/{name}.bam"))
                mapped_reads = int(run_command(f"samtools view -c {bam_path}/{name}_mapped.bam"))

                # Cluster and PCR duplicates (these seem hard-coded in the bash script)
                cluster_duplicates = 5324
                pcr_duplicates = 90000

                # Unique mapped reads (removing duplicates)
                unique_mapping = int(run_command(f"samtools view -c -F 1024 {bam_file}"))

                # Proportions and statistics
                p_trimmed_reads = calculate_proportion(trimmed_reads, total_reads)
                p_mapping = calculate_proportion(mapped_reads, total_reads)
                clonality = calculate_proportion(mapped_reads - unique_mapping, mapped_reads)
                p_cluster_duplicates = calculate_proportion(cluster_duplicates, mapped_reads)
                p_pcr_duplicates = calculate_proportion(pcr_duplicates, mapped_reads)
                endogenous = calculate_proportion(mapped_reads, trimmed_reads)
                p_mapped_unique = calculate_proportion(unique_mapping, trimmed_reads)
                efficiency = calculate_proportion(unique_mapping, total_reads)

                # Read length
                read_length_cmd = f"samtools view -F 1024 {bam_file} | awk '{{sum+=length($10)}} END {{print sum/NR}}'"
                read_length = float(run_command(read_length_cmd))

                x_depth = depth[5]
                y_depth = depth[8]
                m_depth = depth[11]
                aut_depth = depth[14]

                # Write to CSV
                writer.writerow([
                    sample_id, sample_age, library, treatment, 
                    total_reads, trimmed_reads, mapped_reads, 
                    cluster_duplicates, pcr_duplicates, unique_mapping, 
                    p_trimmed_reads, p_mapping, clonality, 
                    p_cluster_duplicates, p_pcr_duplicates, endogenous, 
                    p_mapped_unique, efficiency, read_length, 
                    x_depth, y_depth, m_depth, aut_depth
                ])

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 script.py <project_name>")
        sys.exit(1)
    project_name = sys.argv[1]
    main(project_name)
