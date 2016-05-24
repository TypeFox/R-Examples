Inference <-
function(listComplete,G,M,niter,burnin,threshold=0.5){
  Rres=InferenceR(listComplete[1:niter],G,M,niter,burnin,threshold)
  XiRes=InferenceXi(listComplete[(4*niter+1):(4*niter+3)],niter,burnin)
  ARes=InferenceA(listComplete[(niter+1):(2*niter)],niter,burnin)
  MuRes=InferenceMu(listComplete[(2*niter+1):(3*niter)],niter,burnin)
  SdRes=InferenceSd(listComplete[(3*niter+1):(4*niter)],niter,burnin)
  out=list(R=Rres,Xi=XiRes,A=ARes,Mu=MuRes,Sd=SdRes)
  return(out)
}
