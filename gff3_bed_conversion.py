# Open the input GFF3 file and the output BED file
with open('yeast.gff3', 'r') as gff3_file, open('yeast.bed', 'w') as bed_file:
    # Loop through each line in the GFF3 file
    for line in gff3_file:
        # Skip comment lines that start with '#'
        if line.startswith('#'):
            continue

        # Split the line into columns by tab character
        columns = line.strip().split('\t')

        # Check if there are enough columns to avoid index errors
        if len(columns) < 9:
            continue  # This needs to be inside the loop

        # Check if the feature type (column 3) is a gene
        if columns[2] == 'gene':
            # Get the chromosome name (column 1)
            chrom = columns[0]

            # Get the start and end positions (columns 4 and 5), and convert start to 0-based
            try:
                start = str(int(columns[3]) - 1)
                end = columns[4]
            except ValueError:
                # Handle the case where start/end values are invalid
                continue

            # Set a default name for the gene and look for actual gene name in the attributes (column 9)
            gene_name = 'gene'
            attributes = columns[8].split(';')
            for attribute in attributes:
                if attribute.startswith('Name='):
                    gene_name = attribute.split('=')[1]
                    break

            # Get the strand information (+ or -) (column 7)
            strand = columns[6]

            # Write the BED entry (chromosome, start, end, name, score, strand) to the output file
            bed_file.write(f"{chrom}\t{start}\t{end}\t{gene_name}\t0\t{strand}\n")
