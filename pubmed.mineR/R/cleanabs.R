setGeneric("cleanabs", function(object) standardGeneric("cleanabs"));
setMethod("cleanabs","Abstracts",function(object){temp1 =  which(object@Abstract!="NONE"); temp=new("Abstracts", Journal=object@Journal[temp1], Abstract=object@Abstract[temp1], PMID=object@PMID[temp1]);return(temp)})
