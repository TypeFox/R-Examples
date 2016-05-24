remp<-function(n=NULL,obs=NULL){
  if(is.null(obs)) stop("y is empty.")
  if(class(obs)!="numeric") stop("obs is not a numeric vector.")
  if(n<=0) stop("n must be greater than zero.")
  outpt<-rep(NA,n)
  x<-data.frame(ros=sort(obs[!is.na(obs)]))
  U<-runif(n,0,1)
  x$s<-1:length(x[,1])
  x$GRs<-(x$s-1)/(length(x[,1])-1)
  x$st<-x$s/(length(x[,1])-1)
  for(j in 1:n){
    if(any(x$GRs==U[j])){
      outpt[j]<-x$ros[which(x$GRs==U[j])] 
    } 
    if(all(x$GRs!=U[j])){
      temp<-x[which(U[j]>x$GRs & U[j]<x$st):as.numeric(which(U[j]>x$GRs & U[j]<x$st)+1),]
      outpt[j]<-(length(x[,1])-1)*(temp$ros[2]-temp$ros[1])*(U[j]-temp$GRs[1])+temp$ros[1]
    }
  }
  return(outpt)
}
