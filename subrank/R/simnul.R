simnul <-
function(sampsize,dimension,subsampsizes,sampnum,KL=TRUE,nbsafe=5)
{
  lrs=matrix(ncol=length(subsampsizes),nrow=sampnum)
  scarcities=matrix(ncol=length(subsampsizes),nrow=sampnum)
  for (e in 1:sampnum)
  {
    for (s in 1:length(subsampsizes))
    {
      simdata=rnorm(sampsize*dimension)
      nboot=nbsafe*subsampsizes[s]^dimension
      cop=corc0(simdata,sampsize,dimension,subsampsizes[s],nboot,42)
	  # the coercion to numeric avoids integer overflow, because the same coercion is done implicitally for subsampsizes[s]
      nbootreel=as.numeric(cop[subsampsizes[s]^dimension+2])
      cop=cop[1:subsampsizes[s]^dimension]
      if (KL)
      {
        RV=sum(cop*log(cop),na.rm=TRUE)/(nbootreel*subsampsizes[s])-log(nbootreel*subsampsizes[s])+dimension*log(subsampsizes[s])
      } else
      {
        RV=sum(cop^2)/((nbootreel*subsampsizes[s])^2)-subsampsizes[s]^(-dimension)
      }
      lrs[e,s]=RV
      scarcities[e,s]=sum(cop==0)/(nbootreel*subsampsizes[s])
    }
  }
  return(list(lrs=lrs,scarcities=scarcities))
}
