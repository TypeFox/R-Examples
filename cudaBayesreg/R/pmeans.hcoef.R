# pmeans.hcoef=function(x,burnin=trunc(.1*R),...){
pmeans.hcoef=function(x,burnin=trunc(.1*R))
{
#
# arrays of draws of coefs in hier models
#    3 dimensional arrays:        unit x var x draw
#
    X=x
    if(mode(X) == "list") stop("list entered \n Possible Fixup: extract from list \n")
    if(mode(X) !="numeric") stop("Requires numeric argument \n")
    d=dim(X)
    if(length(d) !=3) stop("Requires 3-dim array \n") 
    nunits=d[1]
    nvar=d[2]
    R=d[3]
    if(R < 100) {cat("fewer than 100 draws submitted \n"); return(invisible())}
    # posterior means for each var 
    # par(las=1)
    pmeans=matrix(0,nrow=nunits,ncol=nvar)
    for(i in 1:nunits) pmeans[i,]=apply(X[i,,(burnin+1):R],1,mean)
    names=as.character(1:nvar)
    attributes(pmeans)$class="hcoef.post"
    for(i in 1:nvar) names[i]=paste("Posterior Means of Coef ",i,sep="")
    return(pmeans)
}

