Crsk <- function(t,ic) {cbind(t,ic)}
Cens <- function(t,ic) {cbind(t,ic)}

crr.wt.KM <- function(t,ic,z,variance,var.conservative) {
  var.name <- colnames(z) 
  
  n <- length(t)
  ncov <- ncol(as.matrix(z))
  tmp <- cbind(t, ic, z)
  sort.tmp <- tmp[order(tmp[,1]),]
  
  ts <- sort.tmp[,1]
  ics <- sort.tmp[,2]
  zs <- matrix(0, nrow=n, ncol=ncov)
  for(i in 1:ncov) zs[,i] <- sort.tmp[,i+2]
  count_event <- length(unique(ts[which(ics==1)]))  	
  
  betaest <- vector(mode="numeric",length=ncov)
  betasd  <- vector(mode="numeric",length=ncov)
  time    <- vector(mode="numeric",length=count_event)
  a10     <- vector(mode="numeric",length=count_event)
  a10sd   <- vector(mode="numeric",length=count_event)
  Wlambda <- vector(mode="numeric",length=1)
  Wbeta   <- vector(mode="numeric",length=1)
  if(variance) {
    Wlambda <- vector(mode="numeric",length=count_event*n)
    Wbeta   <- vector(mode="numeric",length=n*ncov)
  }
  conv    <- 1
  
  out <- .C("FGweight_KM", as.single(ts), as.integer(ics), as.single(t(zs)), 
            as.integer(n), as.integer(ncov), beta=as.single(betaest), 
            beta_sd=as.single(betasd), time=as.single(time), 
            a10=as.single(a10), a10sd=as.single(a10sd), 
            Wlambda=as.single(Wlambda), Wbeta=as.single(Wbeta),
            conv=as.single(conv),as.single(as.numeric(variance)),
            as.single(as.numeric(var.conservative)))
  if(variance) {
    res <- list(weight="KM",
                varname=var.name,
                converge=as.numeric(out$conv),
                beta=out$beta[1:ncov],  
                beta_se=sqrt(out$beta_sd[1:ncov]),
                time=out$time[1:count_event], 
                a10=out$a10[1:count_event], 
                a10se=out$a10sd[1:count_event], 
                W_lambda=out$Wlambda[1:(count_event*n)], 
                W_beta=out$Wbeta[1:(n*ncov)])
  }
  else {
    res <- list(weight="KM",
                varname=var.name,
                converge=as.numeric(out$conv),
                beta=out$beta[1:ncov],  
                beta_se=sqrt(out$beta_sd[1:ncov]),
                time=out$time[1:count_event], 
                a10=out$a10[1:count_event], 
                a10se=out$a10sd[1:count_event], 
                W_lambda=out$Wlambda[1], 
                W_beta=out$Wbeta[1])    
  }
  class(res) <- "crrwt"
  return(res)
}


crr.wt.COX <- function(t,ic,z,zc,variance,var.conservative) {
  var.name <- colnames(z) 
  
  if(is.null(zc)) {
    cat("\n  Warning: Covariate set for censoring is missing. \n  Fit a Cox model for the censoring distribution using terms in argument 'z' as the covariate set\n\n")
    zc <- z
	n <- length(t)
    ncov <- ncol(as.matrix(z))
    ncovc <- ncol(as.matrix(zc))
    tmp <- cbind(t, ic, z, zc)
    sort.tmp <- tmp[order(tmp[,1]),]
    
    ts <- sort.tmp[,1]
    ics <- sort.tmp[,2]
    zs <- matrix(0, nrow=n, ncol=ncov)
    zcs <- matrix(0, nrow=n, ncol=ncovc)
    for(i in 1:ncov) zs[,i] <- sort.tmp[,i+2]
    for(i in 1:ncovc) zcs[,i] <- sort.tmp[,i+2+ncov]
    count_event <- length(unique(ts[which(ics==1)]))  	
    
    betaest <- vector(mode="numeric",length=ncov)
    betasd <- vector(mode="numeric",length=ncov)
    time <- vector(mode="numeric",length=count_event)
    a10 <- vector(mode="numeric",length=count_event)
    a10sd <- vector(mode="numeric",length=count_event)
    Wlambda <- vector(mode="numeric",length=1)
    Wbeta   <- vector(mode="numeric",length=1)
    if(variance) {
      Wlambda <- vector(mode="numeric",length=count_event*n)
      Wbeta <- vector(mode="numeric",length=n*ncov)
    }
    conv <- 1
    censdet <- 1 

    out <- .C("FGweight_COX", as.single(ts), as.integer(ics), as.single(t(zs)), 
              as.single(t(zcs)), as.integer(n), as.integer(ncov), 
              as.integer(ncovc), beta=as.single(betaest), 
              beta_sd=as.single(betasd), time=as.single(time), 
              a10=as.single(a10), a10sd=as.single(a10sd), 
              Wlambda=as.single(Wlambda), Wbeta=as.single(Wbeta),
              conv=as.single(conv),censdet=as.single(censdet),
              as.single(as.numeric(variance)),
              as.single(as.numeric(var.conservative)))
    
    if(as.numeric(out$censdet)==0) {
      cat("Using COX weights had a convergence problem.\nKaplan-Meier weights were used instead.\n\n")
      res <- crr.wt.KM(t,ic,z,variance,var.conservative)      
      return(res)
    }
    
    if(variance) {
      res <- list(weight="COX",
                  varname=var.name,
                  converge=as.numeric(out$conv),
                  cens.det=as.numeric(out$censdet),
                  beta=out$beta[1:ncov],  
                  beta_se=sqrt(out$beta_sd[1:ncov]),
                  time=out$time[1:count_event], 
                  a10=out$a10[1:count_event],
                  a10se=out$a10sd[1:count_event], 
                  W_lambda=out$Wlambda[1:(count_event*n)], 
                  W_beta=out$Wbeta[1:(n*ncov)])
    }
    else {
      res <- list(weight="COX",
                  varname=var.name,
                  converge=as.numeric(out$conv),
                  cens.det=as.numeric(out$censdet),
                  beta=out$beta[1:ncov],  
                  beta_se=sqrt(out$beta_sd[1:ncov]),
                  time=out$time[1:count_event], 
                  a10=out$a10[1:count_event],
                  a10se=out$a10sd[1:count_event], 
                  W_lambda=out$Wlambda[1], 
                  W_beta=out$Wbeta[1])
    }
    class(res) <- "crrwt"
    return(res)
  }
  else { ## when zc is specified
	n <- length(t)
    ncov <- ncol(as.matrix(z))
    ncovc <- ncol(as.matrix(zc))
    tmp <- cbind(t, ic, z, zc)
    sort.tmp <- tmp[order(tmp[,1]),]
    
    ts <- sort.tmp[,1]
    ics <- sort.tmp[,2]
    zs <- matrix(0, nrow=n, ncol=ncov)
    zcs <- matrix(0, nrow=n, ncol=ncovc)
    for(i in 1:ncov) zs[,i] <- sort.tmp[,i+2]
    for(i in 1:ncovc) zcs[,i] <- sort.tmp[,i+2+ncov]
    count_event <- length(unique(ts[which(ics==1)]))	
    
    betaest <- vector(mode="numeric",length=ncov)
    betasd <- vector(mode="numeric",length=ncov)
    time <- vector(mode="numeric",length=count_event)
    a10 <- vector(mode="numeric",length=count_event)
    a10sd <- vector(mode="numeric",length=count_event)
    Wlambda <- vector(mode="numeric",length=1)
    Wbeta   <- vector(mode="numeric",length=1)
    if(variance) {
      Wlambda <- vector(mode="numeric",length=count_event*n)
      Wbeta <- vector(mode="numeric",length=n*ncov)
    }
    conv <- 1
    censdet <- 1
    
    out <- .C("FGweight_COX", as.single(ts), as.integer(ics), 
              as.single(t(zs)), as.single(t(zcs)), as.integer(n), 
              as.integer(ncov), as.integer(ncovc), 
              beta=as.single(betaest), beta_sd=as.single(betasd), 
              time=as.single(time), a10=as.single(a10), 
              a10sd=as.single(a10sd), Wlambda=as.single(Wlambda), 
              Wbeta=as.single(Wbeta), conv=as.single(conv),censdet=as.single(censdet),
              as.single(as.numeric(variance)),
              as.single(as.numeric(var.conservative)))

    if(as.numeric(out$censdet)==0) {
      cat("Using COX weights had a convergence problem.\nKaplan-Meier weights were used instead.\n\n")
      res <- crr.wt.KM(t,ic,z,variance,var.conservative)      
      return(res)
    }
    
    if(variance) {
      res <- list(weight="COX",
                  varname=var.name,
                  converge=as.numeric(out$conv),
                  cens.det=as.numeric(out$censdet),
                  beta=out$beta[1:ncov],  
                  beta_se=sqrt(out$beta_sd[1:ncov]),
                  time=out$time[1:count_event], 
                  a10=out$a10[1:count_event],
                  a10se=out$a10sd[1:count_event], 
                  W_lambda=out$Wlambda[1:(count_event*n)], 
                  W_beta=out$Wbeta[1:(n*ncov)])
    }
    else {
      res <- list(weight="COX",
                  varname=var.name,
                  converge=as.numeric(out$conv),
                  cens.det=as.numeric(out$censdet),
                  beta=out$beta[1:ncov],  
                  beta_se=sqrt(out$beta_sd[1:ncov]),
                  time=out$time[1:count_event], 
                  a10=out$a10[1:count_event],
                  a10se=out$a10sd[1:count_event], 
                  W_lambda=out$Wlambda[1], 
                  W_beta=out$Wbeta[1])      
    }
    
    class(res) <- "crrwt"
    return(res)
  }  
}


crr.wt.KM.str <- function(t,ic,z,strata.var,variance,var.conservative) {
  var.name <- colnames(z) 
  
  # Pre-processing strata.var  
  if(is.null(strata.var)) {
    cat("\n  Warning: Strata variable is missing. \n  Regular Kaplan-Meier weight was used instead.\n\n")
    res <- crr.wt.KM(t,ic,z)
    return(res)
  }
  strata.std <- vector("numeric",length=length(strata.var))
  strata.sort <- sort(unique(strata.var))
  nstrata <- length(unique(strata.var))
  for(i in 1:length(strata.var)) {
    for(j in 1:length(strata.sort)) {
      if(strata.var[i]==strata.sort[j]) {
        strata.std[i] <- j-1
        break
      }
    }
  }
  
  n <- length(t)
  ncov <- ncol(as.matrix(z))
  tmp <- cbind(t, ic, z, strata.std)
  sort.tmp <- tmp[order(tmp[,1]),]
  
  ts <- sort.tmp[,1]
  ics <- sort.tmp[,2]
  zs <- matrix(0, nrow=n, ncol=ncov)
  for(i in 1:ncov) zs[,i] <- sort.tmp[,i+2]
  strata.std <- sort.tmp[,dim(sort.tmp)[2]]
  count_event <- length(unique(ts[which(ics==1)]))    
  
  betaest <- vector(mode="numeric",length=ncov)
  betasd  <- vector(mode="numeric",length=ncov)
  time    <- vector(mode="numeric",length=count_event)
  a10     <- vector(mode="numeric",length=count_event)
  a10sd   <- vector(mode="numeric",length=count_event)
  Wlambda <- vector(mode="numeric",length=1)
  Wbeta   <- vector(mode="numeric",length=1)
  if(variance) {
    Wlambda <- vector(mode="numeric",length=count_event*n)
    Wbeta   <- vector(mode="numeric",length=n*ncov)
  }
  conv    <- 1
  
  out <- .C("FGweight_KM_Strata", as.single(ts), as.integer(ics), as.single(t(zs)), 
            as.integer(n), as.integer(ncov), beta=as.single(betaest), 
            beta_sd=as.single(betasd), time=as.single(time), 
            a10=as.single(a10), a10sd=as.single(a10sd), 
            Wlambda=as.single(Wlambda), Wbeta=as.single(Wbeta),
            conv=as.single(conv),as.integer(nstrata),as.integer(strata.std),
            as.single(as.numeric(variance)),
            as.single(as.numeric(var.conservative)))
  if(variance) {
    res <- list(weight="Stratified KM",
                varname=var.name,
                converge=as.numeric(out$conv),
                beta=out$beta[1:ncov],  
                beta_se=sqrt(out$beta_sd[1:ncov]),
                time=out$time[1:count_event], 
                a10=out$a10[1:count_event], 
                a10se=out$a10sd[1:count_event], 
                W_lambda=out$Wlambda[1:(count_event*n)], 
                W_beta=out$Wbeta[1:(n*ncov)])
  }
  else {
    res <- list(weight="Stratified KM",
                varname=var.name,
                converge=as.numeric(out$conv),
                beta=out$beta[1:ncov],  
                beta_se=sqrt(out$beta_sd[1:ncov]),
                time=out$time[1:count_event], 
                a10=out$a10[1:count_event], 
                a10se=out$a10sd[1:count_event], 
                W_lambda=out$Wlambda[1], 
                W_beta=out$Wbeta[1])    
  }
  class(res) <- "crrwt"
  return(res)  
}


crr.wt.COX.str <- function(t,ic,z,zc,strata.var,variance,var.conservative) {
  var.name <- colnames(z) 
  
  # Pre-processing strata.var   
  if(is.null(strata.var)) {
    cat("\n  Warning: Strata variable is missing. \n  Regular Cox weight was used instead.\n\n")
    res <- crr.wt.COX(t,ic,z,zc)
    return(res)
  }
  strata.std <- vector("numeric",length=length(strata.var))
  strata.sort <- sort(unique(strata.var))
  nstrata <- length(unique(strata.var))
  for(i in 1:length(strata.var)) {
    for(j in 1:length(strata.sort)) {
      if(strata.var[i]==strata.sort[j]) {
        strata.std[i] <- j-1
        break
      }
    }
  }
  
  if(is.null(zc)) {
    cat("\n  Warning: Covariate set for censoring is missing. \n  Fit a stratified Cox model for the censoring distribution using terms in argument 'z' as the covariate set\n\n")
    
	n <- length(t)
	
    # Remove strata variable from regression variables
    strcol <- 0
    for(i in 1:dim(z)[2]) {
      if(sum(z[,i]==strata.var)==n) {
        strcol <- i 
        break
      }
    }
    #      print(i)
    if(strcol > 0) zc <- z[,-strcol]
    else zc <- z
    
    ncov <- ncol(as.matrix(z))
    ncovc <- ncol(as.matrix(zc))
    tmp <- cbind(t, ic, z, zc, strata.std)
    sort.tmp <- tmp[order(tmp[,1]),]
    
    ts <- sort.tmp[,1]
    ics <- sort.tmp[,2]
    zs <- matrix(0, nrow=n, ncol=ncov)
    zcs <- matrix(0, nrow=n, ncol=ncovc)
    for(i in 1:ncov) zs[,i] <- sort.tmp[,i+2]
    for(i in 1:ncovc) zcs[,i] <- sort.tmp[,i+2+ncov]
    strata.std <- sort.tmp[,dim(sort.tmp)[2]]
    count_event <- length(unique(ts[which(ics==1)]))    
    
    betaest <- vector(mode="numeric",length=ncov)
    betasd <- vector(mode="numeric",length=ncov)
    time <- vector(mode="numeric",length=count_event)
    a10 <- vector(mode="numeric",length=count_event)
    a10sd <- vector(mode="numeric",length=count_event)
    Wlambda <- vector(mode="numeric",length=1)
    Wbeta <- vector(mode="numeric",length=1)
    if(variance) {
      Wlambda <- vector(mode="numeric",length=count_event*n)
      Wbeta <- vector(mode="numeric",length=n*ncov)
    }
    conv <- 1
    censdet <- 1 
    
    out <- .C("FGweight_COX_Strata", as.single(ts), as.integer(ics), as.single(t(zs)), 
              as.single(t(zcs)), as.integer(n), as.integer(ncov), 
              as.integer(ncovc), beta=as.single(betaest), 
              beta_sd=as.single(betasd), time=as.single(time), 
              a10=as.single(a10), a10sd=as.single(a10sd), 
              Wlambda=as.single(Wlambda), Wbeta=as.single(Wbeta),
              conv=as.single(conv),censdet=as.single(censdet),
              as.integer(nstrata),as.integer(strata.std),
              as.single(as.numeric(variance)),
              as.single(as.numeric(var.conservative)))
    
    if(as.numeric(out$censdet)==0) {
      cat("Using COX weights had a convergence problem.\n\n")
    }
    
    if(variance) {
      res <- list(weight="Stratified COX",
                  varname=var.name,
                  converge=as.numeric(out$conv),
                  cens.det=as.numeric(out$censdet),
                  beta=out$beta[1:ncov],  
                  beta_se=sqrt(out$beta_sd[1:ncov]),
                  time=out$time[1:count_event], 
                  a10=out$a10[1:count_event],
                  a10se=out$a10sd[1:count_event], 
                  W_lambda=out$Wlambda[1:(count_event*n)], 
                  W_beta=out$Wbeta[1:(n*ncov)])
    }
    else {
      res <- list(weight="Stratified COX",
                  varname=var.name,
                  converge=as.numeric(out$conv),
                  cens.det=as.numeric(out$censdet),
                  beta=out$beta[1:ncov],  
                  beta_se=sqrt(out$beta_sd[1:ncov]),
                  time=out$time[1:count_event], 
                  a10=out$a10[1:count_event],
                  a10se=out$a10sd[1:count_event], 
                  W_lambda=out$Wlambda[1], 
                  W_beta=out$Wbeta[1])     
    }
    
    class(res) <- "crrwt"
    return(res)
  }
  else {
    # Remove strata variable from regression variables
	n <- length(t)
		
    strcol <- 0
    for(i in 1:dim(zc)[2]) {
      if(sum(zc[,i]==strata.var)==n) {
        strcol <- i 
        break
      }
    }
    if(strcol > 0) zc <- zc[,-strcol]
    
    ncov <- ncol(as.matrix(z))
    ncovc <- ncol(as.matrix(zc))
    tmp <- cbind(t, ic, z, zc, strata.std)
    sort.tmp <- tmp[order(tmp[,1]),]
    
    ts <- sort.tmp[,1]
    ics <- sort.tmp[,2]
    zs <- matrix(0, nrow=n, ncol=ncov)
    zcs <- matrix(0, nrow=n, ncol=ncovc)
    for(i in 1:ncov) zs[,i] <- sort.tmp[,i+2]
    for(i in 1:ncovc) zcs[,i] <- sort.tmp[,i+2+ncov]
    strata.std <- sort.tmp[,dim(sort.tmp)[2]]
    count_event <- length(unique(ts[which(ics==1)]))    
    
    betaest <- vector(mode="numeric",length=ncov)
    betasd <- vector(mode="numeric",length=ncov)
    time <- vector(mode="numeric",length=count_event)
    a10 <- vector(mode="numeric",length=count_event)
    a10sd <- vector(mode="numeric",length=count_event)
    Wlambda <- vector(mode="numeric",length=1)
    Wbeta <- vector(mode="numeric",length=1)
    if(variance) {
      Wlambda <- vector(mode="numeric",length=count_event*n)
      Wbeta <- vector(mode="numeric",length=n*ncov)
    }
    conv <- 1
    censdet <- 1 
    
    out <- .C("FGweight_COX_Strata", as.single(ts), as.integer(ics), as.single(t(zs)), 
              as.single(t(zcs)), as.integer(n), as.integer(ncov), 
              as.integer(ncovc), beta=as.single(betaest), 
              beta_sd=as.single(betasd), time=as.single(time), 
              a10=as.single(a10), a10sd=as.single(a10sd), 
              Wlambda=as.single(Wlambda), Wbeta=as.single(Wbeta),
              conv=as.single(conv),censdet=as.single(censdet),
              as.integer(nstrata),as.integer(strata.std),
              as.single(as.numeric(variance)),
              as.single(as.numeric(var.conservative)))
    
    if(as.numeric(out$censdet)==0) {
      cat("Using COX weights had a convergence problem.\n\n")
    }
    
    if(variance) {
      res <- list(weight="Stratified COX",
                  varname=var.name,
                  converge=as.numeric(out$conv),
                  cens.det=as.numeric(out$censdet),
                    beta=out$beta[1:ncov],  
                  beta_se=sqrt(out$beta_sd[1:ncov]),
                  time=out$time[1:count_event], 
                  a10=out$a10[1:count_event],
                  a10se=out$a10sd[1:count_event], 
                  W_lambda=out$Wlambda[1:(count_event*n)], 
                  W_beta=out$Wbeta[1:(n*ncov)])
    }
    else {
      res <- list(weight="Stratified COX",
                  varname=var.name,
                  converge=as.numeric(out$conv),
                  cens.det=as.numeric(out$censdet),
                  beta=out$beta[1:ncov],  
                  beta_se=sqrt(out$beta_sd[1:ncov]),
                  time=out$time[1:count_event], 
                  a10=out$a10[1:count_event],
                  a10se=out$a10sd[1:count_event], 
                  W_lambda=out$Wlambda[1], 
                  W_beta=out$Wbeta[1])      
    }
    
    class(res) <- "crrwt"
    return(res)
  }  
}


crr.wt <- function(formula, data, weight=c("KM","COX","KM.Strata","COX.Strata"), 
                   cens.formula, cause=1, strata.var, variance=TRUE, 
                   var.conservative=FALSE) {
  inobj <- model.frame(formula,data,na.action=NULL)
  n.total <- dim(as.matrix(inobj))[1]
  
  # Pre-processing for formula
  Call <- match.call()
  indx <- match(c("formula","data"),names(Call),nomatch=0)
  if(indx[1]==0) stop("A formula argument is required.")
  temp <- Call[c(1,indx)]
  temp[[1]] <- as.name("model.frame")  
  
  special <- c("strata","cluster","tt")
  if(missing(data)) temp$formula <- terms(formula,special)
  else temp$formula <- terms(formula,special,data=data)

  if(is.R()) m <- eval(temp,parent.frame())
  else m <- eval(temp,sys.parent())
  
  if(nrow(m)==0) stop("No (non-missing) observations")

  Y <- model.extract(m,"response")
  t <- Y[,1]
  ic <- Y[,2]
  n <- length(t)
  n.missing <- n.total - n

  # Factor  - class variables
  fac.indx <- NULL
  for(i in 2:dim(m)[2]) if(is.factor(m[,i])) fac.indx <- c(fac.indx,i)

  fac.len <- NULL
  if(length(fac.indx)>0) {
    fac.label <- NULL
    for(i in 1:length(fac.indx)) {
      tmp.string <- gsub("factor\\(","",colnames(m)[fac.indx[i]])
      tmp.string <- gsub("\\)","",tmp.string)
      tmp.string <- paste(tmp.string," : ",sep="")
      fac.label <- c(fac.label,tmp.string)
    }

    fac.mat <- model.matrix(~factor(m[,fac.indx[1]]))[,-1]
    fac.len <- ncol(as.matrix(fac.mat))
    if(fac.len == 1 && is.vector(fac.mat)) {
      fac.mat <- matrix(fac.mat,ncol=1)
      i <- 1
      colnames(fac.mat) <- paste("factor(m[, fac.indx[",i,"]])",
                                 strsplit(fac.label[1]," ")[[1]][1],sep="")
    }
    if(length(fac.indx)>1) {
      for(i in 2:length(fac.indx)) {
        tmp <- model.matrix(~factor(m[,fac.indx[i]]))[,-1]
        fac.mat <- cbind(fac.mat,tmp)
        fac.len <- c(fac.len,ncol(as.matrix(tmp)))
        if(colnames(fac.mat)[dim(fac.mat)[2]]=="tmp") {
          tmp.string <- paste("factor(m[, fac.indx[",i,"]])",
                              strsplit(fac.label[i]," ")[[1]][1],sep="")
          colnames(fac.mat)[dim(fac.mat)[2]] <- tmp.string
        }
      }
    }

    exp.label <- NULL
    for(i in 1:length(fac.label)) {
      exp.label <- c(exp.label,rep(fac.label[i],fac.len[i]))
    }

    fac.names <- colnames(fac.mat)
    for(i in 1:length(fac.names)) {
      name.split <- strsplit(fac.names[i],")")
      name.common <- paste(name.split[[1]][1],")",sep="")
      name.common <- gsub("\\(","\\\\\\(",name.common)
      name.common <- gsub("\\)","\\\\\\)",name.common)
      name.common <- gsub("\\[","\\\\\\[",name.common)
      name.common <- gsub("\\]","\\\\\\]",name.common)

      colnames(fac.mat)[i] <- gsub(name.common,exp.label[i],fac.names[i])
    }
  
    m <- cbind(m,fac.mat)
    m <- m[,-fac.indx]
  }

  z <- as.matrix(m[,2:dim(m)[2]],nc=dim(m[,2:dim(m)[2]])[2])
  colnames(z) <- colnames(m)[-1]
  var.name <- colnames(z) 
  miss.index <- NULL
  
  # Pre-processing cens
  if(missing(cens.formula)) zc <- NULL
  else {
    if(missing(data)) m2 <- model.frame(cens.formula)
    else m2 <- model.frame(cens.formula,data=data)
    
    ## factor in censoring
    fac.indx <- NULL
    for(i in 2:dim(m2)[2]) if(is.factor(m2[,i])) fac.indx <- c(fac.indx,i)

    if(length(fac.indx) > 0) {
      fac.mat <- model.matrix(~factor(m2[,fac.indx[1]]))[,-1]
      if(is.vector(fac.mat)) fac.mat <- matrix(fac.mat,ncol=1)
      if(length(fac.indx) > 1) {
        for(i in 2:length(fac.indx)) {
          tmp <- model.matrix(~factor(m2[,fac.indx[i]]))[,-1]
          fac.mat <- cbind(fac.mat,tmp)
        }
      }
      m2 <- cbind(m2,fac.mat)
      m2 <- m2[,-fac.indx]
    }
    zc <- as.matrix(m2[,-1])
  }

  
  if(cause > 1) {
    ic <- ifelse(ic==1, 999, ic)  
    ic <- ifelse(ic==cause, 1, ic)
    ic <- ifelse(ic==999, cause, ic)
  }
 
  # Pre-processing strata.var
  z2 <- as.data.frame(z)
  if(missing(strata.var)) strata.var <- NULL
  else if(missing(data)==FALSE && (toString(substitute(strata.var)) %in% names(z2))) {
    strata.var <- eval(parse(text=paste("z2$",
                                        toString(substitute(strata.var)),sep="")))
  }  
  else if(missing(data)==FALSE) {
    strata.var <- eval(parse(text=paste("data$",
                                        toString(substitute(strata.var)),sep="")))
  } 
	
  # Fitting the model
  if(weight=="KM") {
    res <- crr.wt.KM(t,ic,z,variance,var.conservative)
	res$n <- n.total
	res$n.missing <- n.missing
    return(res)
  }
  if(weight=="COX") {
    res <- crr.wt.COX(t,ic,z,zc,variance,var.conservative)
	res$n <- n.total
	res$n.missing <- n.missing
    return(res)
  } 
  if(weight=="KM.Strata") {
    res <- crr.wt.KM.str(t,ic,z,strata.var,variance,var.conservative)
	res$n <- n.total
	res$n.missing <- n.missing
    return(res)
  }  
  if(weight=="COX.Strata") {
    res <- crr.wt.COX.str(t,ic,z,zc,strata.var,variance,var.conservative)
	res$n <- n.total
	res$n.missing <- n.missing
    return(res)
  }
}

predict.crrwt <- function(object,z,...) {
  flag <- 0
  if(sum(object$beta_se)==0) flag <- 1
  
  z2 <- as.matrix(z)
  ncov <- length(object$beta)
  n <- length(object$W_beta)/ncov
  ntime <- length(object$time)
  
  F1 <- matrix(0,nrow=ntime,ncol=dim(z2)[1])
  F1se <- matrix(0,nrow=ntime,ncol=dim(z2)[1])
  
  for(i in 1:dim(z2)[1]) {
    z <- z2[i,]  
    if(length(z) != ncov) {stop("Dimensions do not match!!")}
    
    if(flag==1) {
      ebz <- exp(sum(object$beta*z))
      F1[,i] <- 1-exp(-object$a10*ebz) 
    }
    else {
      tmp <- vector(mode="numeric",length=n)
      sumtmp <- vector(mode="numeric",length=ntime)
      W_F1_i <- matrix(nrow=n,ncol=ntime)
      W_lambda <- matrix(object$W_lambda,nrow=n,byrow=TRUE)
      W_beta <- matrix(object$W_beta,nrow=n,byrow=TRUE)
    
      ebz <- exp(sum(object$beta*z))
      F1[,i] <- 1-exp(-object$a10*ebz)  		
      tmp <- W_beta %*% c(z)
    
      for(j in 1:ntime) {
        W_F1_i[,j] <- ebz*(object$a10[j]*tmp+W_lambda[,j])
      }
    
      W_F1_i2 <- W_F1_i^2
      sumtmp <- colSums(W_F1_i2)
    
      F1se[,i] <- sqrt((1-F1[,i])^2*sumtmp) 
    }
  }

  res <- list(z=z2, time=object$time, F1=F1, F1se=F1se)
  class(res) <- "crrwt.pred"
  return(res)
}

print.crrwt <- function(x, ...) {
  cat("coefficients:\n")
  print(signif(x$beta, 4), ...)
  cat("standard errors:\n")
  print(signif(x$beta_se, 4), ...)
  invisible()
}

plot.crrwt.pred <- function(x, multiple=0, se=0, ...) {
  if(is.null(x)) {
    stop("No object to plot!")
    return
  }
  if(class(x)!="crrwt.pred") {
    stop("Wrong object class!")
  }
  
  nplot <- dim(x$F1)[2]

  if(multiple == 0) {
    par(mfrow=c(1,1))    
    plot(x$time, x$F1[,1],type="s", main="Plot of CIF over time", 
         xlab="Time", ylab="CIF", ylim=c(0,max(x$F1)*1.2))
    if(nplot >= 2) {
      for (i in 2:nplot) {
        points(x$time, x$F1[,i],type="s",lty=i)
      }
    }
    text <- NULL
    for(i in 1:nplot) {
      text <- c(text,paste("Z=(", paste(x$z[i,],collapse=" "), ")", sep=""))
    }  
    legend("topleft","groups",text,cex=0.6,lty=1:nplot,ncol=3)
  }
  else if(multiple == 1) {
    if(nplot==1) par(mfrow=c(1,1))
    if(nplot==2) par(mfrow=c(1,2))
    if(nplot==3) par(mfrow=c(1,3))
    if(nplot==4) par(mfrow=c(2,2))
    if(nplot>4) par(mfrow=c(2,2))
    
    if(se == 0) {
      for (i in 1:nplot) {
        plot(x$time, x$F1[,i],type="s", 
             main=paste("Plot of CIF for\n Z=(", paste(x$z[i,],collapse=" "), ")", sep=""), 
             xlab="Time", ylab="CIF",ylim=c(0,1))
      }
    }
    else {
      for (i in 1:nplot) {
        plot(x$time, x$F1[,i],type="s", 
             main=paste("Plot of CIF for\n Z=(", paste(x$z[i,],collapse=" "), ")", sep=""), 
             xlab="Time", ylab="CIF",ylim=c(0,1))
        F1.lower <- x$F1[,i]-qnorm(0.975)*x$F1se[,i]
        F1.upper <- x$F1[,i]+qnorm(0.975)*x$F1se[,i]
        F1.lower <- ifelse(F1.lower<0,0,F1.lower)
        F1.upper <- ifelse(F1.upper>1,1,F1.upper)       
        points(x$time, F1.lower,type="s",lty=2)
        points(x$time, F1.upper,type="s",lty=2)
      }      
    }
  }
  else cat("No such multiple value!")
}

summary.crrwt <- function(object, ...) {
  nvar <- length(object$beta)
  res <- matrix(0,nrow=nvar,ncol=7)
  colnames(res) <- c("beta","se(beta)","exp(beta)","95% CI lower",
                     "95% CI upper","Z","p-value") 
  
  rownames(res) <- object$varname
  res[,1] <- object$beta
  res[,2] <- object$beta_se
  res[,3] <- exp(object$beta)
  res[,4] <- exp(object$beta-qnorm(.975)*res[,2])
  res[,5] <- exp(object$beta+qnorm(.975)*res[,2])
  res[,6] <- res[,1]/res[,2]
  res[,7] <- round(2*(1-pnorm(abs(res[,6]))),4)
  
  cat("\n Summary results: \n\n")
  cat(" Number of total observations:    ", object$n,"\n")
  cat(" Number of missing observations:  ", object$n.missing,"\n")  
  cat(" Number of observations used:     ", object$n-object$n.missing,"\n\n")
  
  print(round(res,3))
  
  if(object$weight=="KM") {
    cat("\n Censoring distribution was modeled using the Kaplan-Meier estimator \n")
  }
  if(object$weight=="COX") {
    cat("\n Censoring distribution was modeled via Cox regression \n")
  }
  if(object$weight=="Stratified KM") {
    cat("\n Censoring distribution was modeled using the stratified Kaplan-Meier estimator \n")
  }
  if(object$weight=="Stratified COX") {
    cat("\n Censoring distribution was modeled via stratified Cox regression \n")
  }
  invisible()
}
