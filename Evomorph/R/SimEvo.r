SimEvo <- function(vari,consensusvec,resids,ngen,fsamp){
  if (!(is.matrix(resids)))
  {stop("Residuals must be matrix type")}
  if (!(is.vector(consensusvec)))
  {stop("Consensus shape must be vector type")}
  if (vari<=0)
  {stop("Variation coeficient must be nonzero positive integer")}
  if (!(is.vector(consensusvec)))
  {stop("Consensus shape must be in vectorized form")}
  if (missing(ngen))
  {ngen=1000000
  print("Number of generations not provided. Default = 1M years")}
  if (missing(fsamp))
  {fsamp=1000
  print("Frequency of samples not provided. Sample each 1000 years")}
  results=list();
  bh=0;
  bha=0;
  P=cov(resids);
  pca=svd(P);
  w=as.matrix(pca$d);
  u=as.matrix(pca$u);
  for (n in 1:ngen){
  evo=rnorm(ncol(u),mean=0,sd=vari);
  bh=bh+w*evo;
  if(n%%fsamp==0){
  locations=t(bh)%*%solve(u)+consensusvec;
  model=matrix(locations,nrow=(ncol(u)/2),ncol=2,byrow=T);
  results[[n/fsamp]]=model;
  }
  }
  return(results)
}