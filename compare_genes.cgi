#!/usr/bin/python3

# Import necessary modules
import cgi
import os
import re  # Import regex module for parsing

# Tell the web server that this is HTML output
print("Content-Type: text/html\n\n")

# Paths to the files (reference and predicted data)
genbank_file = "/var/www/html/dmurph64/midterm/sequence.gb"
predicted_file = "/var/www/html/dmurph64/midterm/predicted_genes.gbk"


def extract_gene_coordinates(file_path):
    coordinates = []
    with open(file_path, "r") as file:
        in_features = False
        coord_block = ""

        for line in file:
            line = line.strip()

            # Detect the start of the FEATURES section
            if line.startswith("FEATURES"):
                in_features = True

            if in_features:
                if line.startswith("CDS"):
                    coord_block = line  # Start capturing the CDS coordinates
                elif line.startswith("/") or line == "//":  # End of CDS feature block
                    # Adjust to capture both 'complement' and 'join' features
                    match = re.search(r'(complement\()?(\d+)\.\.(\d+)', coord_block)
                    if match:
                        start = int(match.group(2))
                        end = int(match.group(3))
                        strand = '-' if match.group(1) else '+'
                        coordinates.append((start, end, strand))
                    coord_block = ""  # Reset the block after parsing
                else:
                    coord_block += " " + line  # Continue capturing multi-line block

    return coordinates

# Function to compare reference and predicted gene coordinates
def compare_genes(ref_genes, pred_genes):
    summary = {
        "reference_count": len(ref_genes),
        "predicted_count": len(pred_genes),
        "exact_matches": 0,
        "five_prime_agreement": 0,
        "three_prime_agreement": 0,
        "no_overlap": 0
    }

    # Create a set of predicted gene coordinates for efficient comparison
    pred_gene_set = {(start, end): strand for start, end, strand in pred_genes}

    for ref_start, ref_end, ref_strand in ref_genes:
        matched = False

        # Check for exact match
        if (ref_start, ref_end) in pred_gene_set:
            summary["exact_matches"] += 1
            matched = True
        else:
            # Check for 5' agreement
            for (pred_start, pred_end), pred_strand in pred_gene_set.items():
                if ref_start == pred_start and ref_strand == pred_strand:
                    summary["five_prime_agreement"] += 1
                    matched = True
                    break
                elif ref_end == pred_end and ref_strand == pred_strand:
                    summary["three_prime_agreement"] += 1
                    matched = True
                    break

        if not matched:
            summary["no_overlap"] += 1

    return summary

# Extract gene coordinates from both files
reference_genes = extract_gene_coordinates(genbank_file)
#### print(f"Reference Genes: {reference_genes}")
predicted_genes = extract_gene_coordinates(predicted_file)
##### print(f"Predicted Genes: {predicted_genes}")

# Compare the genes
gene_summary = compare_genes(reference_genes, predicted_genes)

# Display the result in HTML
print("<html><head><title>Gene Prediction Comparison</title></head>")
print("<body>")
print("<h1>Gene Prediction Comparison</h1>")
print("<h2>Summary Information</h2>")
print(f"<p><strong>Count of Genes in Reference Annotation:</strong> {gene_summary['reference_count']}</p>")
print(f"<p><strong>Count of Predicted Genes:</strong> {gene_summary['predicted_count']}</p>")
print(f"<p><strong>Count of Exact Matches:</strong> {gene_summary['exact_matches']}</p>")
print(f"<p><strong>Count of 5' Agreement but 3' Disagreement:</strong> {gene_summary['five_prime_agreement']}</p>")
print(f"<p><strong>Count of 3' Agreement but 5' Disagreement:</strong> {gene_summary['three_prime_agreement']}</p>")
print(f"<p><strong>Count of Predicted Genes with No Overlap:</strong> {gene_summary['no_overlap']}</p>")
print("</body></html>")

