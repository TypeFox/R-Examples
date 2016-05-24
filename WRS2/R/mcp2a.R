mcp2a <- function(formula, data, est = "mom", nboot = 599){
  #
  if (missing(data)) {
    mf <- model.frame(formula)
  } else {
    mf <- model.frame(formula, data)
  }
  cl <- match.call()
  
  J <- nlevels(mf[,2])
  K <- nlevels(mf[,3])
  
  grp=NA
  con=0
  alpha=.05
  
  JK <- J * K
  
  est <- get(est)
  
  nfac <- tapply(mf[,1], list(mf[,2],mf[,3]), length, simplify = FALSE)
  nfac1 <- nfac[unique(mf[,2]), unique(mf[,3])]    ## reordering factor levels
  data$row <- unlist(alply(nfac1, 1, sequence), use.names = FALSE)
  
  dataMelt <- melt(data, id = c("row", colnames(mf)[2], colnames(mf)[3]), measured = mf[,1])
  dataWide <- cast(dataMelt, as.formula(paste(colnames(dataMelt)[1], "~", colnames(mf)[2], "+", colnames(mf)[3]))) 
  dataWide$row <- NULL
  x <- dataWide
  
  if(is.matrix(x))
    x <- listm(x)
  if(!is.na(grp)) {
    yy <- x
    for(j in 1:length(grp))
      x[[j]] <- yy[[grp[j]]]
  }
  mvec<-NA
  tempn=0
  for(j in 1:JK) {
    xx <- x[[j]]
    x[[j]] <- xx[!is.na(xx)]
    mvec[j]<-est(x[[j]])
    tempn[j]=length(x[[j]])
  }
  nmax=max(tempn)
  #
  # Create the three contrast matrices
  #
  #             results can be duplicated.
  #  Determine nboot if a value was not specified

  bvec<-matrix(NA,nrow=JK,ncol=nboot)
  #print("Taking bootstrap samples. Please wait.")
  for(j in 1:JK){
    #print(paste("Working on group ",j))
    data<-matrix(sample(x[[j]],size=length(x[[j]])*nboot,replace=TRUE),nrow=nboot)
    bvec[j,]<-apply(data,1,est) # J by nboot matrix, jth row contains
    #                          bootstrapped  estimates for jth group
  }
  outvec<-list()
  temp3<-con2way(J,K)
  for(jj in 1:3){
    con<-temp3[[jj]]
    con<-as.matrix(con)
    ncon<-ncol(con)
    # Determine critical values
    if(alpha==.05){
      dvec<-c(.025,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
      if(ncon > 10){
        avec<-.05/c(11:ncon)
        dvec<-c(dvec,avec)
      }}
    if(alpha==.01){
      dvec<-c(.005,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
      if(ncon > 10){
        avec<-.01/c(11:ncon)
        dvec<-c(dvec,avec)
      }}
    if(nmax>80){
      if(alpha==.05){
        dvec<-c(.05,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
        if(ncon > 10){
          avec<-.05/c(11:ncon)
          dvec<-c(dvec,avec)
        }}
      if(alpha==.01){
        dvec<-c(.01,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
        if(ncon > 10){
          avec<-.01/c(11:ncon)
          dvec<-c(dvec,avec)
        }}
      if(alpha != .05 && alpha != .01){
        dvec<-alpha/c(1:ncon)
      }}
    if(alpha != .05 && alpha != .01)dvec<-alpha/c(1:ncon)
    test<-NA
    bcon<-t(con)%*%bvec #ncon by nboot matrix
    tvec<-t(con)%*%mvec
    for (d in 1:ncon){
      test[d]<-sum(bcon[d,]>0)/nboot
      if(test[d]> .5)test[d]<-1-test[d]
    }
    output<-matrix(0,ncon,6)
    dimnames(output)<-list(NULL,c("con.num","psihat","sig.test","sig.crit","ci.lower","ci.upper"))
    temp2<-order(0-test)
    zvec<-dvec[1:ncon]
    sigvec<-(test[temp2]>=zvec)
    if(sum(sigvec)<ncon){
      dd<-ncon-sum(sigvec) #number that are sig.
      ddd<-sum(sigvec)+1
      zvec[ddd:ncon]<-dvec[ddd]
    }
    output[temp2,4]<-zvec
    icl<-round(dvec[ncon]*nboot)+1
    icu<-nboot-icl-1
    for (ic in 1:ncol(con)){
      output[ic,2]<-tvec[ic,]
      output[ic,1]<-ic
      output[ic,3]<-test[ic]
      temp<-sort(bcon[ic,])
      output[ic,5]<-temp[icl]
      output[ic,6]<-temp[icu]
    }
    outvec[[jj]]<-output
  }
  
  cnamesA <- colnames(mf)[2]
  conA <- temp3[[1]]
  dnamesA <- paste0(cnamesA, 1:ncol(conA))
  cnamesB <- colnames(mf)[3]
  conB <- temp3[[2]]
  dnamesB <- paste0(cnamesB, 1:ncol(conB))
  colnames(conB) <- dnamesB
  dnamesAB <- apply(expand.grid(dnamesA, dnamesB), 1, function(ss) paste(ss[1], ss[2], sep = ":")) 
  conAB <- temp3[[3]]
  contrasts <- as.data.frame(cbind(conA, conB, conAB))
  colnames(contrasts) <- c(dnamesA, dnamesB, dnamesAB)
  rownames(contrasts) <- colnames(dataWide)
  
  alpha.crit <- outvec[[1]][, "sig.crit"]
  outA <- list(psihat = outvec[[1]][, "psihat"], conf.int = outvec[[1]][, c("ci.lower", "ci.upper")], p.value = outvec[[1]][, "sig.test"])
  outB <- list(psihat = outvec[[2]][, "psihat"], conf.int = outvec[[2]][, c("ci.lower", "ci.upper")], p.value = outvec[[2]][, "sig.test"])
  outAB <- list(psihat = outvec[[3]][, "psihat"], conf.int = outvec[[3]][, c("ci.lower", "ci.upper")], p.value = outvec[[3]][, "sig.test"])
  effects <- list(outA, outB, outAB)
  names(effects) <- c(colnames(mf)[2:3], paste(colnames(mf)[2], colnames(mf)[3], sep = ":"))
  
  result <- list(effects = effects, contrasts = contrasts, alpha.crit = alpha.crit, call = cl)
  class(result) <- "mcp"
  result
}
