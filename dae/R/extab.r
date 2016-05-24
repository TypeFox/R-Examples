extab <- function(table, index.factors, order="standard")
# result is a variate that receives values from table variate
{
#
# test arguments
#
	if (mode(index.factors) != "list") stop("Must supply a list")
	nfac <- length(index.factors)
	if (!all(sapply(index.factors, inherits, what="factor"))) 
		stop("All elements of list must be factors or ordereds")
	if (nfac != 1)
    if (var(sapply(index.factors, length)) != 0) stop("All factors must be of the same length")
 	which.ord <- pmatch(casefold(order), c("standard", "yates"), nomatch="")
	if (which.ord == "")	stop("order must be either standard or yates")
# standard order
	if (which.ord == "1") counter <- nfac:1
# Yates order
	else if (which.ord== "2") counter <- 1:nfac
  
#form multiradix counter
	kval <- 1
	radix <- rep(1, length(index.factors[[1]])) 
	for (i in counter)
#reassign factor so unused levels removed
	{ f <- factor(index.factors[[i]])
    nlev <- length(levels(f))
    radix <- ((1:nlev)[f]-1)*kval+radix
    kval <- kval*nlev 
	}
  if (kval != length(table)) 
    stop("Table variate length must be product numbers of levels of factors")

#copy values from table into variate

  expvar <- table[radix]
  return(expvar)
}
