#!/usr/bin/env python3

# Import necessary modules
import cgi  # Handles web form and CGI output (prints to browser)
import os   # Handles file system operations

# Function to read and process the FASTA file
def parse_fasta(file_path):
    # Open the FASTA file to read it
    with open(file_path, 'r') as fasta_file:
        entries = []  # This will hold the final data: header + sequence
        header = None  # To store each header line
        sequence = []  # To store the sequence part of the file

        # Loop through each line in the file
        for line in fasta_file:
            line = line.strip()  # Remove extra spaces or new lines
            if line.startswith('>'):  # If the line starts with '>', it is a header
                if header:  # If there is a previous header, save it before moving on
                    entries.append((header, ''.join(sequence)))  # Store the header + sequence
                header = line[1:]  # Save the header (everything after the '>')
                sequence = []  # Reset the sequence for the next entry
            else:
                sequence.append(line)  # Add the line to the sequence

        # After finishing the loop, add the last entry
        if header:
            entries.append((header, ''.join(sequence)))

        return entries  # Return the processed data (header + sequence)

# Function to break the header into the ID and the rest of the description
def process_header(header):
    parts = header.split(" ", 1)  # Split the header into two parts by the first space
    seq_id = parts[0]  # The ID part is the first part
    rest_of_header = parts[1] if len(parts) > 1 else ""  # The rest of the header
    return seq_id, rest_of_header

# Set the file path for the FASTA file
fasta_file = '/var/www/html/dmurph64/unit04/e_coli_k12_dh10b.faa'

# Call the function to parse the file
fasta_entries = parse_fasta(fasta_file)

# Limit the number of entries to display (only show the first 20)
max_entries = 20

# Print the necessary CGI header (for the browser to understand it's HTML)
print("Content-type: text/html\n")

# Start of the HTML content
print("<html>")
print("<head><title>FASTA File Information</title></head>")
print("<body>")
print("<h1>FASTA File Information</h1>")

# Create a table to display the information in a structured format
print("<table border='1'>")
print("<tr><th>ID</th><th>Length</th><th>Rest of Header</th></tr>")

# Loop through the first 20 entries and print them in the table
for i, (header, sequence) in enumerate(fasta_entries[:max_entries]):
    # Process the header to extract ID and the rest
    seq_id, rest_of_header = process_header(header)
    # Get the sequence length
    seq_length = len(sequence)
    # Print each row of the table
    print(f"<tr><td>{seq_id}</td><td>{seq_length}</td><td>{rest_of_header}</td></tr>")

# Close the table and HTML body
print("</table>")
print("</body></html>")
