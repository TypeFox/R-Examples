intertest <- function(x,y,k=5){
  ## x is control, y is treatment
  ## return S, E, V, T, P
  ##Based on Orban and Wolfe 1982 JASA 77: 666-672
  ##Changed their expression (2.2) to be a sum from 0 to m
  m<-length(x)
  n<-length(y)
  phi<-function(h){choose(m*h,k-1)/(n*choose(m,k-1))}
  U<-as.vector(unlist(apply(outer(x,y,"<="),2,sum)))/m
  const<-(n*choose(m,k-1))
  Snm<-sum(phi(U))*const
  ph<-phi((0:m)/m)*const
  phb<-sum(ph)/(m+1)
  ESnm<-n*phb
  varSnm<-(sum(ph*ph)-((m+1)*(phb*phb)))*(n*(n+m+1)/((m+1)*(m+2)))
  out<-matrix(NA,1,5)
  colnames(out)<-c("Score","Exp","Var","Dev","p.value")
  rownames(out)<-"Result"
  out[1,1:3]<-c(Snm,ESnm,varSnm)
  out[1,4]<-(Snm-ESnm)/sqrt(varSnm)
  out[1,5]<-1-pnorm(out[1,4])
  ## out<-as.data.frame(out)
  as.numeric(out)
}
