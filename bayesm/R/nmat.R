nmat=function(vec) 
{
#
#    function to take var-cov matrix in vector form and create correlation matrix 
#       and store in vector form
#
p=as.integer(sqrt(length(vec)))
sigma=matrix(vec,ncol=p)
nsig=1/sqrt(diag(sigma))
return(as.vector(nsig*(t(nsig*sigma))))
}

