coeffvar <-
function(x){
  n<-length(x)
  cv<- sd(x)/mean(x)
  cvstar<- (1+1/(4*n))*cv # bias corrected
  list(cv=cv,"bccv"=cvstar) 
}
