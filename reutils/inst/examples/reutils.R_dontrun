#
# combine esearch and efetch
#
# Download PubMed records that are indexed in MeSH for both 'Chlamydia' and 
# 'genome' and were published in 2013.
query <- "Chlamydia[mesh] and genome[mesh] and 2013[pdat]"

# Upload the PMIDs for this search to the History server
pmids <- esearch(query, "pubmed", usehistory = TRUE)
pmids

\dontrun{
# Fetch the records
articles <- efetch(pmids)

# Use XPath expressions with the #xmlValue() or #xmlAttr() methods to directly
# extract specific data from the XML records stored in the 'efetch' object.
titles <- articles$xmlValue("//ArticleTitle")
abstracts <- articles$xmlValue("//AbstractText")

#
# combine epost with esummary/efetch
#
# Download protein records corresponding to a list of GI numbers.
uid <- c("194680922", "50978626", "28558982", "9507199", "6678417")

# post the GI numbers to the Entrez history server
p <- epost(uid, "protein")

# retrieve docsums with esummary
docsum <- content(esummary(p, version = "1.0"), "parsed")
docsum

# download FASTAs as 'text' with efetch
prot <- efetch(p, retmode = "text", rettype = "fasta")
prot

# retrieve the content from the efetch object
fasta <- content(prot)
}
