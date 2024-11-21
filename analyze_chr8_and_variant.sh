#!/bin/bash

# Ensure tabix indexing for the original VCF file
echo "Indexing centenarian.vcf.gz..."
tabix -p vcf centenarian.vcf.gz

# Extract chromosome 8 data into a new compressed VCF
echo "Extracting chr8 data from centenarian.vcf.gz..."
bcftools view -r chr8 -Oz -o chr8.vcf.gz centenarian.vcf.gz

# Reindex chr8.vcf.gz to ensure efficient querying
echo "Indexing chr8.vcf.gz..."
tabix -p vcf chr8.vcf.gz

# Inspect the data without headers
echo "Inspecting chr8.vcf.gz without headers..."
bcftools view chr8.vcf.gz --no-header | head

# Look for a specific position (1178735) in chr8.vcf.gz
echo "Extracting data for chr8:1178735..."
bcftools view -r chr8:1178735 chr8.vcf.gz

# Extract rs6601606 variant by ID into a separate file
echo "Extracting rs6601606 to rs6601606.vcf..."
bcftools view -i 'ID="rs6601606"' chr8.vcf.gz > rs6601606.vcf

# Extract allele frequency information into a plain text file
echo "Extracting allele frequencies to rs6601606_af.txt..."
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\n' rs6601606_af.vcf > rs6601606_af.txt

# Verify allele frequency for position 11780735
echo "Verifying allele frequencies for chr8:11780735..."
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\n' chr8.vcf.gz | grep "11780735"

# Include additional annotation such as total alleles (AN)
echo "Extracting allele frequencies with total allele count (AN)..."
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%AN\n' rs6601606_af.vcf > rs6601606_af_with_an.txt

# End of script
echo "Analysis complete!"
