weighted.coeffvar <-
function(x,w){
  n<-length(x)
  cv<- sqrt(cov.wt(as.matrix(x,ncol=1),wt=w)$cov)/weighted.mean(x,w)
  cvstar<- (1+1/(4*n))*cv # bias corrected
  list(cv=cv,"bccv"=cvstar) 
  print("Warning: weighting is not properly accounted for in the sample adjustment of bccv!!")
}
