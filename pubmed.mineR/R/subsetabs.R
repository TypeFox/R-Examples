setGeneric("subsetabs", function(object,indices) standardGeneric("subsetabs"));
setMethod("subsetabs","Abstracts",function(object,indices){temp=new("Abstracts", Journal=object@Journal[indices], Abstract=object@Abstract[indices], PMID=object@PMID[indices]);return(temp)})
