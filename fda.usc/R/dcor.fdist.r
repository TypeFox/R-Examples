################################################################################
# Description
# Computes distance correlation for functional data.
#
# Arguments
# D1: class(D1)=matrix: distances of first sample,  class(D1)=fdata: first fdata data sample
# D2: class(D2)=matrix: distances of second sample, class(D2)=fdata: second fdata data sample
# Returns the sample distance correlation
cor.fdist=function(D1,D2=NULL,...){
  if (is.null(D2)) {
        if (is.fdata(D1)) D1=metric.dist(D1,...)
        m1row=rowMeans(D1)#apply(D1,1,mean)
        m1col=colMeans(D1)#apply(D1,2,mean)
        n=nrow(D1)
        ones=rep(1,n)
        pD1=D1-outer(ones,m1row)-outer(m1col,ones)+mean(D1)
        res=sqrt(sum(pD1*pD1))/n
        out<-res/sqrt(res*res) 
  }
  else {
        if (is.fdata(D1)) D1=metric.lp(D1,...)      
        if (is.fdata(D2)) D2=metric.lp(D2,...)                    
        m1row=rowMeans(D1) 
        m2row=rowMeans(D2)
        m1col=colMeans(D1)
        m2col=colMeans(D2)
        n=nrow(D1)
        ones=rep(1,n)
        pD1=D1-outer(ones,m1row)-outer(m1col,ones)+mean(D1)
        pD2=D2-outer(ones,m2row)-outer(m2col,ones)+mean(D2)
        res=sqrt(sum(pD1*pD2))/n
        res1=sqrt(sum(pD1*pD1))/n        
        res2=sqrt(sum(pD2*pD2))/n        
        out<-res/sqrt(res1*res2)
}
return(out)
}
################################################################################
################################################################################


