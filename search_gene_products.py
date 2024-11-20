#!/usr/bin/env python3
import cgi
import mysql.connector  # Ensure the mysql-connector-python library is installed
import html

# Database connection setup 
conn = mysql.connector.connect(
    host='localhost',
    user='dmurph64',
    password='P@ssw0rd',
    database='chado'
)

# Create a cursor object
cursor = conn.cursor()

# Get the search term from the form
form = cgi.FieldStorage()
search_term = form.getfirst("search_term", "").strip()

# Output HTTP header for the HTML response
print("Content-type: text/html\n")

# Start HTML output
print("<html>")
print("<head><title>Gene Product Search Results</title></head>")
print("<body>")
print("<h1>Search Results for '{0}'</h1>".format(html.escape(search_term)))

# SQL query to search for gene products securely using a parameterized query
query = """
    SELECT 
        f.uniquename AS feature, 
        f.feature_id, 
        fp.value AS gene_product_name
    FROM 
        feature f
    JOIN 
        featureprop fp ON f.feature_id = fp.feature_id
    JOIN 
        cvterm cvt ON fp.type_id = cvt.cvterm_id
    WHERE 
        cvt.name = 'gene_product_name'
        AND fp.value LIKE %s
"""

try:
    # Execute query using the user-provided search term (securely preventing SQL injection)
    cursor.execute(query, (f"%{search_term}%",))

    # Fetch the results
    results = cursor.fetchall()

    if results:
        # Display the results in a table
        print("<table border='1'>")
        print("<tr><th>Feature</th><th>Unique Name</th><th>Gene Product</th></tr>")
        for row in results:
            feature, feature_id, gene_product_name = row
            print(
                f"<tr><td>{html.escape(str(feature_id))}</td><td>{html.escape(feature)}</td><td>{html.escape(gene_product_name)}</td></tr>")
        print("</table>")
    else:
        print("<p>No results found.</p>")
except Exception as e:
    print(f"<p>Error occurred: {html.escape(str(e))}</p>")

# Clean up cursor and connection
finally:
    cursor.close()
    conn.close()

# Close HTML output
print("</body></html>")
