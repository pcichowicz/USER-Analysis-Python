import sys

def main(file_path):
    # Initialize variables
    x_reads = y_reads = m_reads = aut_reads = 0
    x_sites = y_sites = m_sites = aut_sites = 0
    x_cov = y_cov = m_cov = aut_cov = 0

    try:
        # Open and read input file
        with open(file_path, 'r') as file:
            for line in file:
                # Split the line into columns
                fields = line.strip().split()
                if len(fields) < 3:
                    continue  # Skip lines that don't have at least 3 columns
                chr_, pos, cov = fields[0], int(fields[1]), int(fields[2])

                # Process by chromosome
                if chr_ == "chrX":
                    x_reads += cov
                    x_sites += 1
                    if cov > 0:
                        x_cov += 1
                elif chr_ == "chrY":
                    y_reads += cov
                    y_sites += 1
                    if cov > 0:
                        y_cov += 1
                elif chr_ == "chrM":
                    m_reads += cov
                    m_sites += 1
                    if cov > 0:
                        m_cov += 1
                else:
                    aut_reads += cov
                    aut_sites += 1
                    if cov > 0:
                        aut_cov += 1
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

    # Print the results
    print(f"{'Chr':<8}{'Sites':<8}{'Depth':<8}")
    print(f"xDepth  {x_cov:<8}{x_reads / x_sites if x_sites > 0 else 0:.2f}")
    print(f"yDepth  {y_cov:<8}{y_reads / y_sites if y_sites > 0 else 0:.2f}")
    print(f"mDepth  {m_cov:<8}{m_reads / m_sites if m_sites > 0 else 0:.2f}")
    print(f"autDepth {aut_cov:<8}{aut_reads / aut_sites if aut_sites > 0 else 0:.2f}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 script.py <input_file>")
        sys.exit(1)
    input_file = sys.argv[1]
    main(input_file)
