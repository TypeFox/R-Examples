wbs.lsw <-
function(y,C_i=tau.fun(y),scales=NULL,M=0,cstar=0.75,lambda=.75) {
    n=length(y)
    J=floor(log(n,2))
    if(is.null(scales)){
      scales<-J-1:floor(3*lambda*log(log(n)))
      max.scale<-min(scales, J-floor(J/2)+1)
    } else{
      scales<-sort(scales, decreasing=T)
      max.scale<-min(scales)
    } 
    if(length(scales)==1) stop(".........Choose at least two scales.........")
    epp<-c()
    for (j in 1:length(scales)) {
      ep = round(max(2*n/2^scales[j], ceiling(sqrt(n)/2)))
      epp = c(epp,ep)
    }
    epp=round(epp/2) 
    z=ews.trans(y,scales=scales)
    dis<-c((n-n/2^(min(scales))+1):n)
    z=z[-dis,]
    u<-uh.wbs(z, scale=scales,epp=epp,del=floor(log(length(y))^2/3),C_i=C_i,M=M,cstar=cstar)
    cp.out=u$breakpoints
    OUT=post.processing(z,del=floor(log(length(y))^2/3),br=cp.out,C_i=C_i,epp=epp,scales=scales)
    suppressWarnings(if (is.na(OUT)) OUT=NULL)
    list(cp.bef=cp.out,cp.aft=OUT)
  }
