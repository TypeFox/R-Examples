read.design<-function(m,Terms,model="aalen"){ ## {{{ 
  mt <- attr(m, "terms")
  intercept <- attr(mt, "intercept")

  XZ<-model.matrix(Terms,m)[, drop = FALSE]

  clusterTerms<- grep("^cluster[(][A-z0-9._:]*",colnames(XZ),perl=TRUE)

  if(length(clusterTerms) == 1){
    l.cols<-length(attributes(XZ)$assign) - 1
    clusters <- as.vector(XZ[,clusterTerms]);
    cols<-attributes(XZ)$assign
    cols <- cols[-clusterTerms]
  } else if(length(clusterTerms) == 0){
    l.cols<-length(attributes(XZ)$assign)
    cols<-attributes(XZ)$assign
    clusters <- NULL;
  } else {
    stop("More than one cluster(-) term was included in the formula")
  }

  if (model=="aalen" || model=="comprisk") semicov <- grep("^const[(][A-z0-9._:]*",colnames(XZ),perl=TRUE)
  if (model=="prop.excess") semicov <- grep("^cox[(][A-z0-9._:]*",colnames(XZ),perl=TRUE)
  if (model=="cox.aalen") semicov <- grep("^prop[(][A-z0-9._:]*",colnames(XZ),perl=TRUE)
  
   ZtermsXZ<-semicov-1

  if (length(semicov)) {
    npar<-FALSE; Zterms<-c();
    for (i in ZtermsXZ) Zterms<-c(Zterms,(1:l.cols)[cols==i]); 
  } else {
    npar<-TRUE;
  }
  Zterms<-semicov # ts 25-6-2008

  if (length(semicov)>0) {
   covnamesX <- colnames(XZ)[-c(Zterms,clusterTerms)];
   px <- length(covnamesX); 
   X<-matrix(XZ[,-c(Zterms,clusterTerms)],ncol=px);
   if (px==length(colnames(X))) colnames(X)<-covnamesX;
   covnamesZ <- colnames(XZ)[Zterms];
   pz <- length(covnamesZ); 
   Z<-matrix(XZ[,Zterms],ncol=pz);
   if (pz==length(colnames(Z))) colnames(Z)<-covnamesX;
  } else if(length(clusterTerms)>0) {
    covnamesX <- colnames(XZ)[-clusterTerms]
    px <- length(covnamesX); 
    X<-matrix(XZ[,-clusterTerms],ncol=px);
    if (px==length(colnames(X))) colnames(X)<-covnamesX;
  } else {
    covnamesX <- colnames(XZ)
    px <- length(covnamesX); 
    X<-matrix(XZ,ncol=px);
    if (px==length(colnames(X))) colnames(X)<-covnamesX;
  }

  px <- ncol(X)
  if (npar == FALSE) pz <- ncol(Z) else pz <- 0
  pxz <- px + pz

  X <- data.matrix(X)
  if (npar == FALSE){
    Z <- data.matrix(Z)
  }else {
    Z <- NULL; covnamesZ<-NULL
  }
  if(length(clusterTerms) > 0){
    clusters <- as.vector(clusters)
  }

  ud<-list(X=X,Z=Z,px=px,pz=pz,npar=npar,
           covnamesX=covnamesX,covnamesZ=covnamesZ,
           clusters=clusters)
   return(ud)
} ## }}} 

read.surv<-function(m,id,npar,clusters,start.time,max.time,model="aalen",silent=0){

  if (attr(m[, 1], "type") == "right") {
    time2 <- m[, 1][, "time"]; time <- rep(0, length(time2))
    status <- m[, 1][, "status"]
  }
  else if (attr(m[, 1], "type") == "counting") {
###    if ((is.null(id)==TRUE) && (silent==0)) 
###      cat("For counting process data, with multiple records the id variable must be set to get correct robust standard errors and simulations\n");
    time <- m[, 1][, 1];time2 <- m[, 1][, 2];status <- m[, 1][, 3];
  }
  else { stop("only right-censored or counting processes data") }

  if (sum(duplicated(time2[status==1]))>0) {
    # cat("Non unique survival times: break ties ! \n")
    # cat("Break ties yourself\n");
    ties<-TRUE; 
    dtimes<-time2[status==1]
    index<-(1:length(time2))[status==1]
    ties<-duplicated(dtimes); nties<-sum(ties); index<-index[ties]
    dt<-diff(sort(time2)); dt<-min(dt[dt>0]); 
    time2[index]<-time2[index]+runif(nties,0,min(0.001,dt/2));
  } else ties<-FALSE; 

  if ((model=="aalen") &  (npar==FALSE)) times<-unique(time2) else times<-time2[status==1];

  times <- c(start.time, times[times>start.time]); times <- sort(times)
  if (is.null(max.time) == TRUE) maxtimes <- max(times) else maxtimes <- max.time
  times<-times[times<=maxtimes];
  if ((npar==FALSE) & (model!="cox.aalen")) times<-c(times,maxtimes)  
  times<-unique(times); Ntimes <- length(times); 

  if (is.null(id) == TRUE) {antpers <- length(time); id <- 0:(antpers - 1);}
  else {pers <- unique(id); antpers <- length(pers);
        id <- as.integer(factor(id, labels = 1:(antpers))) - 1; 
  }

  if (is.null(clusters)== TRUE) {clusters<-id; antclust<-antpers;} else {
    clus<-unique(clusters); antclust<-length(clus); 
    clusters <- as.integer(factor(clusters, labels = 1:(antclust))) - 1;
  }

  ud2<-list(status=status,start=time,stop=time2,antpers=antpers,antclust=antclust,
            times=times,id.call=id,clusters=clusters,cluster)
  return(ud2)
}

check.missing<-function(X,Z,time,time2,status,npar)
{
XZ<-cbind(X,Z)
  if  ((prod(is.na(XZ) == FALSE) == 0) || (prod(is.na(time) == 
    FALSE) == 0) || (prod(is.na(time2) == FALSE) == 0) || 
   (prod(is.na(status) == FALSE) == 0)) cat("Missing values\n")

}

rm.missing<-function(X,Z,time,time2,status,npar,na.rm=TRUE)
{
  if (is.null(Z)==TRUE) XZ<-cbind(X,Z) else XZ<-X; 

  if ((prod(is.na(XZ) == FALSE) == 0) || (prod(is.na(time) == 
             FALSE) == 0) || (prod(is.na(time2) == FALSE) == 0) || 
      (prod(is.na(status) == FALSE) == 0)) {
    if (na.rm) {
      notmissing <- c(rep(1, times = dimcovar[1]))
      i <- 1
      for (i in 1:dimcovar[1]) {
        notmissing[i] <- prod(is.na(XZ[i, ]) == FALSE) }
      notmissing <- notmissing * (is.na(time) == FALSE) * 
        (is.na(time2) == FALSE) * (is.na(status) == FALSE)
      cat(dimcovar[1] - sum(notmissing), " observations ignored because of missing values\n")
      dimcovar[1] <- sum(notmissing)
      XZ<- XZ[notmissing == 1, ]
      X <- X[notmissing == 1, ]
      time <- time[notmissing == 1, ]
      time2 <- time2[notmissing == 1, ]
      if (npar== FALSE) Z <- Z[notmissing == 1, ]
      trisk <- trisk[notmissing == 1]
      status <- status[notmissing == 1]
    }
    else { stop("Missing not allowed unless na.rm=TRUE\n\n") }
  }

  return(list(X=X,Z=Z,time=time,time2=time2,status=status))
}

