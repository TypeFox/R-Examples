"fcab" <-
function(last,nints,yam1,h,x,stdv){
  f <- last*dnorm(h*c(0:nints)+yam1,mean=matrix(rep(x,nints+1),nints+1,length(x),byrow=TRUE),sd=stdv)
  area <- 0.5*h*(2*colSums(f)-f[1,]-f[nrow(f),]) # This is "trap"
  return(area)
}

