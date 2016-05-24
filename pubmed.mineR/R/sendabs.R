setGeneric("sendabs", function(object,x) standardGeneric("sendabs"));
setMethod("sendabs", "Abstracts", function(object,x){write.table(cbind(object@Journal,object@Abstract,object@PMID), file = x,quote=FALSE, sep = "\t", row.names = FALSE, col.names = c("Journal", "Abstract", "PMID")  )})
