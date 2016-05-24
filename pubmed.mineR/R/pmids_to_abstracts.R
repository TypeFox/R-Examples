pmids_to_abstracts = function(x,abs){ids=NULL; for ( i in 1:length(x)) ids = c(ids,which(abs@PMID == x[i])); return(subsetabs(abs,ids))}
