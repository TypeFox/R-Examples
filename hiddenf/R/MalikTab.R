MalikTab <-
function(r,c,N=1000){
  rows=as.factor(rep(1:r,each=c))
  cols<-as.factor(rep(seq(1:c),r))
  Tcsim<-c()
  for(i in 1:N){
    ysim<-rnorm((r*c),0,1)
    modsim<-lm(ysim~rows + cols)
    rsim<-resid(modsim)
    rmatsim<-matrix(rsim,nrow=length(rsim),ncol=1)
    kmeansim<-kmeans(x=rmatsim,centers=3,nstart=100)
    assnsim<-kmeansim$cluster
    modclussim<-lm(ysim~rows + cols +as.factor(assnsim))
    amodclussim<-anova(modclussim)
    Tcsim[i]<- (amodclussim[3,2]/amodclussim[3,1])/(amodclussim[4,2]/amodclussim[4,1])}
    return(list(Tcsim=sort(Tcsim),q=c(r=r,c=c,quantile(Tcsim,c(.99,.95,.9)))))
  }
