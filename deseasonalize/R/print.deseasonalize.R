`print.deseasonalize` <-
function(x, ...){
dspar <- x$dspar
if (nrow(dspar)==1) out<-dspar else out <- dspar[rownames(dspar)=="*",]
names(out) <- colnames(dspar)
print.default(out)
}

