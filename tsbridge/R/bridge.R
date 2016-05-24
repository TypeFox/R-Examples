bridge <-
function(w1,w2,r0,tol,verbose=FALSE){
  if(sum(colnames(w1) %in% c("q1","q2","l"))!=3 | sum(colnames(w2) %in% c("q1","q2","l"))!=3)
    stop("columns of w1 and w2 must be named q1, q2 and l")
  
  rhat.fn<-function(l1,l2,r){
    l1.star<-max(l1)
    n1<-length(l1)
    n2<-length(l2)
    s1<-n1/(n1+n2)
    s2<-n2/(n1+n2)
    list(r=log(n1/n2) + log(sum( 1 / (s1+s2*exp(r-l2) ))) +
           l1.star - log(sum( 1 / (s1*exp(l1-l1.star)+s2*exp(r-l1.star)) )), n1=n1, n2=n2 )
    #		list(r=log(n1/n2) + log(sum( exp(l2-l2.star) / (s1*exp(l2-l2.star)+s2*exp(r-l2.star)))) +
    #			  l1.star - log(sum( 1 / (s1*exp(l1-l1.star)+s2*exp(r-l1.star)) )), n1=n1, n2=n2 )
  }
  old.r<-0
  r<-r0
  while(abs(old.r-r)>tol){
    old.r<-r
    rhat<-rhat.fn(w1$l, w2$l, old.r)
    r<-rhat$r
    if(verbose==TRUE)	print(r)
  }
  cat(c("sample l1 =" ,rhat$n1, "\n"))
  cat(c("sample l2 =" ,rhat$n2, "\n"))
  return(r)
}
