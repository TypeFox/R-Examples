summary.marqLevAlg <- function(object,digits=8,...){
x <- object
if (!inherits(x, "marqLevAlg")) stop("use only with \"marqLevAlg\" objects")

cat(" \n")
cat("Values:", "\n")
cat(" \n")
id <- 1:length(x$b)
indice <- rep(id*(id+1)/2)
se <-sqrt(x$v[indice])
wald <- (x$b/se)**2
z <- abs(qnorm((1 + .95)/2))
tmp <- data.frame("coef"=format(round(x$b,digits)),"SE coef"=format(round(se,digits)),"Wald"=format(wald,digits),"P-value"=format.pval(1 - pchisq(wald, 1),digits=digits,eps=0.0001))

print(tmp,row.names=F)

}

