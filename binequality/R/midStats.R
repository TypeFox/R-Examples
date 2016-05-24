midStats <-
function(data){
  means<-c()
  medians<-c()
  ginis<-c()
  theils<-c()
  cvs<-c()
  IDs<-c()
  MLDs <-c ()
  for(i in unique(data$ID)){
    use.i <- which(data$ID == i)
    dat.i <- rep(data$mids[use.i], times = data$hb[use.i])
    mu.i<-mean(dat.i,na.rm=TRUE)
    med.i<-median(dat.i, na.rm=TRUE)
    gini.i<-Gini(dat.i)
    theil.i<-Theil(dat.i)
    cv.i<- var.coeff(dat.i,square=FALSE)
    MLD.i <- MLD(na.omit(dat.i))
    
    IDs<-c(IDs,i)
    means<-c(means,mu.i)
    medians<-c(medians,med.i)
    ginis<-c(ginis,gini.i)
    theils<-c(theils,theil.i)
    cvs<-c(cvs,cv.i)
    MLDs <-c (MLDs, MLD.i)
  }
  out<-data.frame(IDs,means,medians,ginis,theils,cvs, MLDs)
  colnames(out)<-c('ID','mean','median','gini','theil','cv','MLD')
  return(out)
}
