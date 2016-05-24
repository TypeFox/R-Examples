mcp2atm<-function(formula, data, tr = 0.2){
  
  if (missing(data)) {
    mf <- model.frame(formula)
  } else {
    mf <- model.frame(formula, data)
  }
  cl <- match.call()
  
  J <- nlevels(mf[,2])
  K <- nlevels(mf[,3])
  
  alpha=.05
  grp=NA
  op=F
  JK <- J * K
  
  
  nfac <- tapply(mf[,1], list(mf[,2],mf[,3]), length, simplify = FALSE)
  nfac1 <- nfac[unique(mf[,2]), unique(mf[,3])]    ## reordering factor levels
  data$row <- unlist(alply(nfac1, 1, sequence), use.names = FALSE)
  
  dataMelt <- melt(data, id = c("row", colnames(mf)[2], colnames(mf)[3]), measured = mf[,1])
  dataWide <- cast(dataMelt, as.formula(paste(colnames(dataMelt)[1], "~", colnames(mf)[2], "+", colnames(mf)[3]))) 
  dataWide$row <- NULL
  x <- dataWide
  
  if(is.matrix(x))
    x <- listm(x)
  if(!is.na(grp[1])) {
    yy <- x
    x<-list()
    for(j in 1:length(grp))
      x[[j]] <- yy[[grp[j]]]
  }
  for(j in 1:JK) {
    xx <- x[[j]]
    x[[j]] <- xx[!is.na(xx)] # Remove missing values
  }
  #
  
  for(j in 1:JK){
    temp<-x[[j]]
    temp<-temp[!is.na(temp)] # Remove missing values.
    x[[j]]<-temp
  }
  # Create the three contrast matrices
  temp<-con2way(J,K)
  conA<-temp$conA
  conB<-temp$conB
  conAB<-temp$conAB
  if(!op){
    Factor.A<-lincon1(x,con=conA,tr=tr,alpha=alpha)
    Factor.B<-lincon1(x,con=conB,tr=tr,alpha=alpha)
    Factor.AB<-lincon1(x,con=conAB,tr=tr,alpha=alpha)
  }
  All.Tests<-NA
  if(op){
    Factor.A<-NA
    Factor.B<-NA
    Factor.AB<-NA
    con<-cbind(conA,conB,conAB)
    All.Tests<-lincon1(x,con=con,tr=tr,alpha=alpha)
  }
  
  cnamesA <- colnames(mf)[2]
  dnamesA <- paste0(cnamesA, 1:ncol(conA))
  cnamesB <- colnames(mf)[3]
  dnamesB <- paste0(cnamesB, 1:ncol(conB))
  colnames(conB) <- dnamesB
  dnamesAB <- apply(expand.grid(dnamesA, dnamesB), 1, function(ss) paste(ss[1], ss[2], sep = ":")) 
  contrasts <- as.data.frame(cbind(conA, conB, conAB))
  colnames(contrasts) <- c(dnamesA, dnamesB, dnamesAB)
  rownames(contrasts) <- colnames(dataWide)
  
  outA <- list(psihat = Factor.A[[3]][, "psihat"], conf.int = Factor.A[[3]][, c("ci.lower", "ci.upper")], p.value = Factor.A[[3]][, "p.value"])
  outB <- list(psihat =  Factor.B[[3]][, "psihat"], conf.int = Factor.B[[3]][, c("ci.lower", "ci.upper")], p.value = Factor.B[[3]][, "p.value"])
  outAB <- list(psihat = Factor.AB[[3]][, "psihat"], conf.int = Factor.AB[[3]][, c("ci.lower", "ci.upper")], p.value = Factor.AB[[3]][, "p.value"])
  effects <- list(outA, outB, outAB)
  names(effects) <- c(colnames(mf)[2:3], paste(colnames(mf)[2], colnames(mf)[3], sep = ":"))
  
  result <- list(effects = effects, contrasts = contrasts, call = cl)
  class(result) <- "mcp"
  result
}
