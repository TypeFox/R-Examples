LLRA <- function(X, W, mpoints, groups, baseline=NULL, itmgrps=NULL,...)
{
  if(missing(mpoints)) stop("Please specify the number of time points. If there are none, you might want to try PCM() or LPCM().")
  Xprep <- llra.datprep(X,mpoints,groups,baseline)
  itmgrps <- rep(1:Xprep$nitems) 
  groupvec <- Xprep$assign.vec
  pplgrps <- length(Xprep$grp_n)
  if(missing(W)) W <- build_W(Xprep$X,length(unique(itmgrps)),mpoints,Xprep$grp_n,groupvec,itmgrps)
  fit <- LPCM(Xprep$X,W,mpoints=mpoints,groupvec=groupvec,sum0=FALSE)
  refg <- unique(names(which(groupvec==max(groupvec))))
  out <- c(fit,"itms"=Xprep$nitems,"refGroup"=refg)
  out$call <- match.call()
  class(out) <- c("llra","Rm","eRm")
  cat("Reference group: ",refg,"\n\n")
  return(out)
}
