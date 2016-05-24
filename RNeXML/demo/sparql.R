library(rrdf)
library(phytools)
library(RNeXML)


nexml <- nexml_read(system.file("examples/primates.xml", package="RNeXML"))

# Extract the RDF graph from the nexml
rdf <- get_rdf(system.file("examples/primates.xml", package="RNeXML"))
saveXML(rdf, "rdf_meta.xml") # rrdf requires a file name, so we must write the XML out first
graph <- load.rdf("rdf_meta.xml")


# fetch the NCBI URI for the taxon that has rank 'Order', i.e. the root of the primates. The dot operator
# '.' between clauses implies a join, in this case
root <- sparql.rdf(graph, 
    "SELECT ?uri WHERE { 
        ?id <http://rs.tdwg.org/ontology/voc/TaxonConcept#rank> <http://rs.tdwg.org/ontology/voc/TaxonRank#Order> . 
        ?id <http://rs.tdwg.org/ontology/voc/TaxonConcept#toTaxon> ?uri    
    }"               
)

# Define a function to get the name 
get_name <- function(id) {
    max <- length(nexml@otus[[1]]@otu)
    for(i in 1:max) {
        if ( nexml@otus[[1]]@otu[[i]]@id == id ) {
            label <- nexml@otus[[1]]@otu[[i]]@label
            label <- gsub(" ","_",label)
            return(label)
        }
    }
}

# define a recursive function to build newick
recurse <- function(node){
  
    # fetch the taxonomic rank and id string
    rank_query <- paste0(
        "SELECT ?rank ?id WHERE {
            ?id <http://rs.tdwg.org/ontology/voc/TaxonConcept#toTaxon> <",node,"> .
            ?id <http://rs.tdwg.org/ontology/voc/TaxonConcept#rank> ?rank
          }")
    result <- sparql.rdf(graph, rank_query)
    
    # get the local ID, strip URI part
    id <- result[2]
    id <- gsub("^.+#", "", id, perl = TRUE)
    
    # if rank is terminal, return the name
    if (result[1] == "http://rs.tdwg.org/ontology/voc/TaxonRank#Species") {
        return(get_name(id))
    }
    
    # recurse deeper
    else {
        child_query <- paste0(
            "SELECT ?uri WHERE {
                ?id <http://www.w3.org/2000/01/rdf-schema#subClassOf> <",node,"> .
                ?id <http://rs.tdwg.org/ontology/voc/TaxonConcept#toTaxon> ?uri
            }")
        children <- sparql.rdf(graph, child_query)
        
        # the newick can be made to contain interior node labels by inserting get_name(id), before the sep="" argument
        return(paste("(", 
                     paste(sapply(children, recurse), 
                           sep = ",", collapse = "," ), 
                     ")",  sep = "", collapse = ""))
    }
}

# build the tree and visualize it
newick <- paste(recurse(root), ";", sep = "", collapse = "")
tree <- read.newick(text = newick)
collapsed <- collapse.singles(tree)
plot(collapsed, type = "cladogram")
