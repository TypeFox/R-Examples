pd <-
function(d,iv=1,jv=2) {
  m<-length(d[1,])
  if(m<2){
    print("data set must have more than 1 variable (column).....")
    return(NULL)
  }
  names<-names(d)
  if(length(names)<2){
    print("data set must be a data.frame with named variables.......")
    return(NULL)
  }
  V1<-d[,iv]
  V2<-d[,jv]
  d<-data.frame(V1,V2)
  d <- d[complete.cases(d),]
  aScore<-ma(d)
  n<-length(d[,1])
  rd <- rwt(d)       
  logLi <- c(0.0)
  kvec <- c(aScore$marginalKW,aScore$jointKW)
  weight <- c(aScore$weight)
  dist <- rep(0.0,times=n*n)
  ret <- .C("getDistribution", np=as.integer(n), rdp=as.integer(rd),   
            kvecp=as.double(kvec), weightp=as.double(weight), distp=as.double(dist))
  optDist <- ret$distp
  dim(optDist)<-c(n,n)
  contour(optDist,
          main=paste("Density Contour Plot, n = ",toString(n),
                     "\n A = ",toString(round(aScore$A,2)),
                     " rawA = ",toString(round(aScore$rawA,2)),
                     "\n jointKW = ",toString(round(aScore$jointKW,2)),
                     " marginalKW = ",toString(round(aScore$marginalKW,2)),
                     " weight = ",toString(round(aScore$weight,2))
                     ),
          xlab=names[iv],
          ylab=names[jv]
          )
  return(dist)
}
