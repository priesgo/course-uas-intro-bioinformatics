


## Practical use case 1: gathering data from a genomic database

We are going to fetch the assembled sequence of SARS-CoV-2 isolated in Wuhan in 2019. This will work as our reference genome for SARS-CoV-2.

### Step 1: fetch the desired accession

Open the REST API interface Swagger for the ENA browser https://www.ebi.ac.uk/ena/browser/api/swagger-ui.html

Beware, there are two REST APIs for ENA.

Perform a **GET search** with the following parameters parameters:
- result=sequence, ie: we want assembled sequences
- query=tax_id=2697049 AND collection_date=2019-12-01
	- The tax id 2697049 corresponds to the SARS-CoV-2 virus. It is a bit daunting but the API works with these tax identifiers instead of scientific names of organisms.
	- We want to get one of the first sequences extracted at Wuhan hence we fetch the assemblies from the 1st December 2019
- format=json, ie: JSON is a computer friendly format very easy to work with in Python

Execute the query!

```
[ 
	{ 
		"description": "Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome.", 
		"collection_date": "Dec-2019", 
		"accession": "MN908947", 
		"tax_id": "2697049" 
	} 
]
```

Query URL
https://www.ebi.ac.uk/ena/portal/api/search?result=sequence&query=tax_id%3D2697049%20AND%20collection_date%3D2019-12-01&limit=10&format=json

So, we got our accession **MN908947**! Now we need to fetch the actual sequence for this accession. Such sequences are stored in FASTA format.

## Step 2: fetch the sequence

To fetch the sequence we need to use the other REST API (I am pretty sure these two will become one eventually).

https://www.ebi.ac.uk/ena/browser/api/swagger-ui.html

Perform a **POST fasta** with the accession as only parameter.

Execute the query!

```
>ENA|MN908947|MN908947.3 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome.
ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT
GTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACT
CACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATC
[...]
```

Query URL (beware this being a POST it won't work to copy this in your browser as the browser does GET requests by default).
https://www.ebi.ac.uk/ena/browser/api/fasta?accessions=MN908947

## Step 3: now bring it all to Python

Perform the same queries programatically:
```
import requests
import json

# prepare the query
tax_id = 2697049
collection_date = "2019-12-01"
query = "tax_id={} AND collection_date={}".format(tax_id, collection_date)
url = "https://www.ebi.ac.uk/ena/portal/api/search?result=sequence&query={}&limit=10&format=json".format(query)

# execute the query
response = requests.get(url)
print(response.content)

# now parse the response and extract the accession
response_json = json.loads(response_content)

# get the first element from the response, ie: there is only one but it is in a list
# and extract the value of the field "accession"
accession = response_json[0].get("accession")
print(accession)

# now let's query the other REST API
url2 = "https://www.ebi.ac.uk/ena/browser/api/fasta?accessions={}".format(accession)

# execute the second query
response2 = requests.post(url2)

# print the first five lines of the response
print(str(response2.content).split("\\n")[0:5])

# let's write the sequence into a FASTA file in mode w=write and b=binary
fasta = open("reference_genome.fasta", mode="wb")
fasta.write(response2.content)
fasta.close()
```

### Exercises

- How many bases does the SARS-CoV-2 sequence we gathered have? 
- 29903
```
sum([len(x) for x in str(response2.content).split("\\n")[1:-1]])
```

- If the GC content is calculated as the ratio of bases being either G or C, calculate the GC content of the SARS-CoV-2 genome.
- 37.97 %
```
g_count = sum([1 if y=='G' else 0 for x in str(response2.content).split("\\n")[1:-1] for y in x])
c_count = sum([1 if y=='C' else 0 for x in str(response2.content).split("\\n")[1:-1] for y in x])
float(g_count + c_count) / 29903
```




