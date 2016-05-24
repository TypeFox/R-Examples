setGeneric("subabs", function(object,start, end) standardGeneric("subabs"));
setMethod("subabs","Abstracts",function(object,start, end){temp=new("Abstracts", Journal=object@Journal[start:end], Abstract=object@Abstract[start:end], PMID=object@PMID[start:end]);return(temp)})
