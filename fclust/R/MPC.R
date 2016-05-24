MPC <-
function (U)
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
k=ncol(U)
if (k==1) 
{ 
cat("There is only k=1 cluster: the MPC index is not computed",fill=TRUE)
mod.part.coeff=NA
}
else 
mod.part.coeff=1-k/(k-1)*(1-PC(U))
return(mod.part.coeff)
}
