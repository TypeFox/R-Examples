mlgarch <-
function(y, arch=1, garch=1, xreg=NULL,
  initial.values=NULL, lower=NULL, upper=NULL,
  nlminb.control=list(), vcov=TRUE, objective.penalty=NULL,
  solve.tol=.Machine$double.eps, c.code=TRUE)
{
  #check/change arguments:
  if(is.null(arch)){ arch <- 0 }
  if(is.null(garch)){ garch <- 0 }
  if(arch < garch) stop("garch order cannot be greater than arch order, since estimation is via the varma representation")
  if(arch > 1) stop("Sorry, arch order cannot be greater than 1 in the current version of mlgarch")

  #zoo:
  y <- as.zoo(y)
  y <- na.trim(y)
  y.index <- index(y)
  y <- coredata(y)
  y.colnames <- colnames(y)
  colnames(y) <- NULL

  #xreg:
  if(!is.null(xreg)){
    if(NROW(xreg)!=NROW(y)) stop("NROW(xreg) must equal NROW(y)")
    xreg <- as.matrix(coredata(as.zoo(xreg)))
    xreg.colnames <- colnames(xreg)
    colnames(xreg) <- NULL
  }

  #rows and dimensions:
  aux <- list()
  aux$y <- y; #rm(y) in the future?
  aux$y.index <- y.index
  aux$n <- NROW(y)
  aux$m <- NCOL(y)
  if(!is.null(xreg)){
    aux$xreg.m <- NCOL(aux$xreg)
    aux$xreg <- xreg; #rm(xreg) in the future?
  }

  #orders:
  aux$maxpq <- max(arch,garch)
  aux$nmaxpq <- aux$n + aux$maxpq
  aux$ar <- aux$maxpq
  aux$ma <- garch

  #zeros and more:
  aux$yiszero <- matrix(as.numeric(y==0),aux$n,aux$m)
  aux$yiszeroadj <- rbind(matrix(0,aux$maxpq,aux$m),aux$yiszero)
  aux$yanyrowiszero <- as.numeric(rowSums(aux$yiszero)>0)
  aux$yanyrowiszeroadj <- as.numeric(rowSums(aux$yiszeroadj)>0)
  aux$yzerowhichrows <- which(aux$yanyrowiszero==1)
  aux$yanyrowiszeron <- sum(aux$yanyrowiszero)
  aux$ynonzerorowsn <- aux$n - aux$yanyrowiszeron
  aux$yzeron <- colSums(aux$yiszero)
  aux$ynonzeron <- aux$n - aux$yzeron
  aux$yzerowhere <- list()
  for(i in 1:aux$m){
    if(aux$yzeron[i] > 0){
      aux$yzerowhere[[i]] <- which(y[,i] == 0)
    }else{ aux$yzerowhere[[i]] <- numeric(0) }
  } #end aux$zerowhere list
  y2 <- y^2
  miny2 <- rep(NA,aux$m)
  for(i in 1:aux$m){
    if(aux$yzeron[i] > 0){
      miny2[i] <- min(y2[-aux$yzerowhere[[i]],i])
      y2[aux$yzerowhere[[i]],i] <- miny2[i]
    }
  }
  aux$lny2 <- log(y2)
  aux$Elny2 <- colMeans(aux$lny2)
  for(i in 1:aux$m){
    if(aux$yzeron[i] > 0){
      aux$Elny2[i] <- mean(aux$lny2[-aux$yzerowhere[[i]],i])
    }
  }
  aux$lny2adj <- aux$lny2
  if(aux$maxpq > 0){
    aux$lny2adj <- rbind(t(matrix(rep(aux$Elny2,aux$maxpq),aux$m,aux$maxpq)),
      aux$lny2adj)
#    for(i in 1:aux$maxpq){
#      aux$lny2adj <- rbind(aux$Elny2,aux$lny2adj)
#    }
  }
  rm(y2) #remove y2 for memory-efficiency reasons

  #indices:
  aux$const.indx <- 1:aux$m
  if(aux$ar>0){
    aux$ar.indx <- c(aux$m+1):c(aux$m+aux$ar*aux$m^2)
  }else{
    aux$ar.indx <- 0
  }
  if(aux$ma>0){
    aux$ma.indx <- c(max(aux$ar.indx)+1):c(max(aux$ar.indx)+aux$ma*aux$m^2)
  }else{
    aux$ma.indx <- 0
  }
  if(is.null(aux$xreg)){
    aux$xreg.k <- 0
    aux$xreg.indx <- 0
  }else{
    aux$xreg.k <- aux$xreg.m * aux$m
    aux$xreg.indx <- c(max(aux$m,aux$ar.indx,aux$ma.indx)+1):c(max(aux$m,aux$ar.indx,aux$ma.indx)+aux$xreg.k)
  }
  aux$sigma2u.indx <- c(max(aux$m,aux$ar.indx,aux$ma.indx,aux$xreg.indx)+1):c(max(aux$m,aux$ar.indx,aux$ma.indx,aux$xreg.indx)+aux$m)
  if(aux$m>1){
    aux$cov.indx <- c(max(aux$sigma2u.indx)+1):c(max(aux$sigma2u.indx)+(aux$m^2-aux$m)/2)
  }else{
    aux$cov.k <- 0
    aux$cov.indx <- 0
  }

  #initial values:
  if(is.null(initial.values)){
    if(aux$ma > 0){
      ma.initvals <- as.vector( diag(rep(-0.8/aux$ma, aux$m)) )
    }else{ ma.initvals <- NULL }
    if(aux$ar > 0){
      ar.initvals <- as.vector( diag(rep(0.9/aux$ar, aux$m)) )
    }else{ ar.initvals <- NULL }
    #future: check ar.initvals for stability?
    if(is.null(aux$xreg)){
      xreg.initvals <- numeric(0)
      aux$xregMeans <- 0
    }else{
      xreg.initvals <- rep(0.01, aux$xreg.k)
      aux$xregMeans <- as.numeric(colMeans(aux$xreg))
    }
    IminPhi1 <- diag(rep(1,aux$m))
    if(aux$ar > 0){
      IminPhi1 <- IminPhi1 - matrix(ar.initvals,aux$m,aux$m)
    }
    const.initvals <- IminPhi1%*%cbind(aux$Elny2)
    if(aux$xreg.k > 0){
      const.initvals <- const.initvals - matrix(xreg.initvals,aux$m,aux$xreg.m) %*% aux$xregMeans
    }
    const.initvals <- as.vector(const.initvals)
    sigma2u.initvals <- rep(4.94,aux$m)
    if(aux$m>1){
      cov.initvals <- var(aux$lny2)[lower.tri(var(aux$lny2))]
    }
    initial.values <- c(const.initvals, ar.initvals, ma.initvals,
      xreg.initvals, sigma2u.initvals, cov.initvals)
  }else{
    if( length(initial.values)!=max( aux$cov.indx ) ){
      stop("length(initial.values) not equal to no. of parameters to be estimated")
    } #end check length(initial.values)
  } #end if..else..
  aux$initial.values <- initial.values

  #upper bounds:
  if(is.null(upper)){
    upper <- rep(Inf,aux$m) #constant
    upper <- c(upper, rep(1-.Machine$double.eps,aux$ar*aux$m^2)) #ar parameters
    upper <- c(upper, rep(1-.Machine$double.eps,aux$ma*aux$m^2)) #ma parameters
    upper <- c(upper, rep(Inf,aux$xreg.k))
    upper <- c(upper, rep(Inf,aux$m)) #variances
    upper <- c(upper, rep(Inf,(aux$m^2-aux$m)/2)) #covariances
  }else{
    if( length(upper)!=length(initial.values) )
      stop("length(upper) not equal to length(initial.values)")
  }
  aux$upper <- upper

  #lower bounds:
  if(is.null(lower)){
    lower <- rep(-Inf,aux$m) #constant
    lower <- c(lower, rep(-1+.Machine$double.eps,aux$ar*aux$m^2)) #ar parameters
    lower <- c(lower, rep(-1+.Machine$double.eps,aux$ma*aux$m^2)) #ma parameters
    lower <- c(lower, rep(-Inf,aux$xreg.k))
    lower <- c(lower,rep(0,aux$m)) #variances
    lower <- c(lower, rep(-Inf,(aux$m^2-aux$m)/2)) #covariances
  }else{
    if( length(lower)!=length(initial.values) )
      stop("length(lower) not equal to length(initial.values)")
  }
  aux$lower <- lower

  #misc:
  aux$c.code <- c.code
  aux$solve.tol <- solve.tol
  aux$verboseRecursion <- FALSE
#  aux$yzeroadj <- c(rep(1,max(1,aux$maxpq)), aux$yzero)
#  aux$zerosaux <- rep(0,max(aux$nmaxpq, aux$n+1))
  if(is.null(objective.penalty)){
    aux$objective.penalty <- mlgarchObjective(initial.values, aux)
  }else{
    aux$objective.penalty <- objective.penalty
  }

  #estimate:
  objective.f <- function(pars, x=aux){ -mlgarchObjective(pars,x) }
  est <- nlminb(initial.values, objective.f, lower=lower,
    upper=upper, control=nlminb.control)
  est$objective <- -est$objective
  names(est)[2] <- "objective.varma"

  #parameters:
  uadj <- mlgarchRecursion1(as.numeric(est$par), aux)
  Elnz2 <- rep(NA,aux$m)
  for(i in 1:aux$m){
    if(aux$yzeron[i] > 0){
      uadjtmp <- uadj[-aux$yzerowhere[[i]],i]
    }else{
      uadjtmp <- uadj[,i]
    }
    Elnz2[i] <- -log(mean(exp(uadjtmp - mean(uadjtmp))))
  }
  rm(uadjtmp) #remove object
  parMlgarch <- Elnz2
  namesMlgarch <- paste("Elnz2no",1:aux$m,sep="")
  parVarma <- est$par
  matIndices <- matrix(NA,aux$m,aux$m)
  for(j in 1:aux$m){
    for(i in 1:aux$m){
      matIndices[i,j] <- paste(i,j,sep="")
    }
  }
  namesVarma <- matIndices[lower.tri(matIndices)]
  namesVarma <- c(diag(matIndices),namesVarma)
  namesVarma <- paste("cov",namesVarma,sep="")
#  parMlgarch <- c(parVarma[aux$cov.indx],parMlgarch)
#  parMlgarch <- c(parVarma[aux$sigma2u.indx],parMlgarch)
#  namesMlgarch <- c(namesVarma,namesMlgarch)
  if(aux$xreg.k > 0){
    xregIndices <- matrix(NA,aux$m,aux$xreg.m)
    for(j in 1:aux$xreg.m){
      for(i in 1:aux$m){
#        xregIndices[i,j] <- paste(i,j,sep="")
        xregIndices[i,j] <- paste(i,"no",j,sep="")
      }
    }
    xregNames <- paste("xreg",as.character(xregIndices),sep="")
    namesVarma <- c(xregNames, namesVarma)
    namesMlgarch <- c(xregNames, namesMlgarch)
    parMlgarch <- c(est$par[aux$xreg.indx], parMlgarch)
  }
  if(aux$ma > 0){
    namesVarma <- c(paste("ma",as.character(matIndices),
      ".1",sep=""), namesVarma)
    namesMlgarch <- c(paste("garch",as.character(matIndices),
      ".1",sep=""), namesMlgarch)
    parMlgarch <- c(-parVarma[aux$ma.indx], parMlgarch)
  }
  if(aux$ar > 0){
    namesVarma <- c(paste("ar",as.character(matIndices),
      ".1",sep=""), namesVarma)
    namesMlgarch <- c(paste("arch",as.character(matIndices),
      ".1",sep=""), namesMlgarch)
    if(aux$ma > 0){
      parTmp <- c(parVarma[aux$ar.indx]+parVarma[aux$ma.indx])
    }else{
      parTmp <- parVarma[aux$ar.indx]
    }
    parMlgarch <- c(parTmp, parMlgarch)
  }
  namesVarma <- c(paste("intercept",1:aux$m,sep=""),namesVarma)
  namesMlgarch <- c(paste("intercept",1:aux$m,sep=""),namesMlgarch)
  if(aux$ma > 0){
    tmpMlgarch <- parVarma[aux$const.indx] - (diag(rep(1,aux$m))+matrix(parVarma[aux$ma.indx],aux$m,aux$m))%*%Elnz2
  }else{
    tmpMlgarch <- parVarma[aux$const.indx] - diag(rep(1,aux$m))%*%Elnz2
  }
  parMlgarch <- c(tmpMlgarch, parMlgarch)
  names(parVarma) <- namesVarma
  names(parMlgarch) <- namesMlgarch
  est$par <- parMlgarch
  est <- c(list(date=date(), par.varma=parVarma), est)

  #vcov matrix:
  if(vcov){
    hessian.varma <- -optimHess(as.numeric(parVarma), objective.f)
    colnames(hessian.varma) <- namesVarma
    rownames(hessian.varma) <- namesVarma
    vcov.varma <- solve(-hessian.varma, tol=aux$solve.tol)
    est <- c(list(aux=aux, hessian.varma=hessian.varma,
      vcov.varma=vcov.varma, vcov.mlgarch=NULL), est)
    est$vcov.mlgarch <- vcov.mlgarch(est)
  }

  #out:
  if(is.null(est$aux)){
    est <- c(list(aux=aux),est)
  }
  class(est) <- "mlgarch"
  return(est)
}
