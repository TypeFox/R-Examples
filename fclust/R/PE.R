PE <-
function (U, b)
{
if (missing(U))
stop("The membership degree matrix U must be given")
if (is.null(U))
stop("The membership degree matrix U is empty")
U=as.matrix(U)
if (any(is.na(U)))
stop("The membership degree matrix U must not contain NA values")
if (!is.numeric(U)) 
stop("The membership degree matrix U is not numeric")
if (missing(b))
{
b=exp(1)
}
if (!is.numeric(b)) 
{
b=exp(1)
cat("The logarithmic base b is not numeric: the default value b=exp(1) will be used ",fill=TRUE)
}
if (b<=1) 
{
b=exp(1)
cat("The logarithmic base b must be >1: the default value b=exp(1) will be used ",fill=TRUE)
}
U[U<.Machine$double.eps]=.Machine$double.eps 
part.ent=-sum(U*log(U,b))/nrow(U)
return(part.ent)
}
