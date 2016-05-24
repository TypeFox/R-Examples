lgarch <-
function(y, arch=1, garch=1, xreg=NULL,
  initial.values=NULL, lower=NULL, upper=NULL,
  nlminb.control=list(), vcov=TRUE, method=c("ls","ml","cex2"),
  mean.correction=FALSE, objective.penalty=NULL,
  solve.tol=.Machine$double.eps, c.code=TRUE)
{
  #check/change arguments:
  if(arch < garch) stop("garch order cannot be greater than arch order, since estimation is via the arma representation")
  if(arch > 1) stop("Sorry, arch order cannot be greater than 1 in the current version of lgarch")
  method <- match.arg(method)
  if(method=="cex2"){ mean.correction <- TRUE }
  if(mean.correction==TRUE && method=="ls" && arch==0 && garch==0 && is.null(xreg) )
    stop("This combination is not possible. Try setting mean.correction=FALSE")

  #zoo and xts specific:
  y.name <- deparse(substitute(y))
  y <- as.zoo(cbind(y))
  if( is.null(y.name)){ y.name <- colnames(y)[1] }
  if( y.name[1] =="" ){ y.name <- "y" }
  y <- na.trim(y)
  y.n <- NROW(y)
  y.index <- index(y)
  t1 <- y.index[1]
  t2 <- y.index[y.n]
  y <- coredata(y)
  if(NCOL(y) > 1){
    stop("Dependent variable not 1-dimensional")
  }else{
    y <- y[,1]
  }

  ##xreg:
  if(!is.null(xreg)){
    xreg <- as.zoo(cbind(xreg))
    xreg.names <- colnames(xreg)
    if(is.null(xreg.names)){
      xreg.names <- paste("xreg", 1:NCOL(xreg), sep="")
    }
    if(any(xreg.names == "")){
      missing.colnames <- which(xreg.names == "")
      for(i in 1:length(missing.colnames)){
        xreg.names[i] <- paste("xreg", i, sep="")
      }
    }
    xreg <- window(xreg, start=t1, end=t2)
    xreg <- cbind(coredata(xreg))
  }

  ##begin aux list:
  aux <- list()
  aux$method <- method
  aux$y <- y
  aux$y.index <- y.index
  aux$n <- length(y)
  aux$yzero <- as.numeric(y != 0)
  aux$ynonzeron <- sum(aux$yzero)
  aux$yzeron <- aux$n - aux$ynonzeron
  if(aux$yzeron > 0) aux$yzerowhere <- which(y == 0)
  y2 <- y^2
  if(aux$yzeron > 0){
    miny2 <- min(y2[-aux$yzerowhere])
    y2[aux$yzerowhere] <- miny2
    #delete?:
    #lny2 <- log(y2)
    aux$lny2 <- log(y2)
    aux$Elny2 <- mean(aux$lny2[-aux$yzerowhere])
  }else{
    aux$lny2 <- log(y2)
    aux$Elny2 <- mean(aux$lny2)
  }
  if(mean.correction){ aux$lny2mc <- aux$lny2 - aux$Elny2 }

  #orders:
  aux$maxpq <- max(arch,garch)
  aux$nmaxpq <- aux$n + aux$maxpq
  aux$ar <- aux$maxpq
  aux$ma <- garch
  if(aux$ar>0){ aux$ar.indx <- 2:c(aux$ar+1) }else{ aux$ar.indx <- 0 }
  if(aux$ma>0){
    aux$ma.indx <- c(max(aux$ar.indx)+1):c(max(aux$ar.indx)+aux$ma)
  }else{
    aux$ma.indx <- 0
  }
  if(is.null(xreg)){
    aux$xreg.k <- 0
  }else{
    aux$xreg <- xreg
    aux$xreg.k <- ncol(aux$xreg)
    aux$xreg.indx <- c(max(1,aux$ar.indx,aux$ma.indx)+1):c(max(1,aux$ar.indx,aux$ma.indx)+aux$xreg.k)
  }
  if(method!="ls"){
    aux$sigma2u.indx <- 1 + aux$ar + aux$ma + aux$xreg.k + 1
  }

  #initial values:
  if(is.null(initial.values)){
    if(is.null(xreg)){
      xreg.initvals <- numeric(0)
      xregMean <- 0
    }else{
      xreg.initvals <- rep(0.01, aux$xreg.k)
      xregMean <- mean(aux$xreg %*% xreg.initvals)
    }
    ma.initvals <- rep(-0.8/aux$ma, aux$ma)
    ar.initvals <- rep(0.9/aux$ar, aux$ar)
    #to do: check ar.initvals for stability?
    constant <- (1-sum(ar.initvals))*aux$Elny2 - xregMean
    sigma2u <- NULL
    if(method=="ml"){ sigma2u <- 4.94 }
    if(method=="cex2"){ sigma2u <- -1.27 }
    initial.values <- c(constant, ar.initvals, ma.initvals,
      xreg.initvals, sigma2u)
  }else{
    max.indx <- max( c(1,aux$ar.indx,aux$ma.indx,aux$xreg.indx,aux$sigma2u.indx) )
    if( length(initial.values)!=max.indx ){
      stop("length(initial.values) not equal to no. of parameters to be estimated")
    } #end check length
  } #end if..else..

  #upper bounds:
  if(is.null(upper)){
    upper <- c(Inf, rep(1-.Machine$double.eps,aux$ar),
      rep(1-.Machine$double.eps,aux$ma), rep(Inf, aux$xreg.k))
    if(method!="ls"){ upper <- c(upper,Inf) }
  }else{
    if( length(upper)!=length(initial.values) )
      stop("length(upper) not equal to length(initial.values)")
  }
  aux$upper <- upper

  #lower bounds:
  if(is.null(lower)){
    lower <- c(-Inf, rep(.Machine$double.eps-1,aux$ar),
      rep(.Machine$double.eps-1,aux$ma), rep(-Inf, aux$xreg.k))
    if(method=="ml"){ lower <- c(lower,0) }
    if(method=="cex2"){ lower <- c(lower,-Inf) }
  }else{
    if( length(lower)!=length(initial.values) )
      stop("length(lower) not equal to length(initial.values)")
  }
  aux$lower <- lower

  #misc:
  aux$c.code <- c.code
  aux$solve.tol <- solve.tol
  aux$mean.correction <- mean.correction
  aux$verboseRecursion <- FALSE
  aux$yzeroadj <- c(rep(1,max(1,aux$maxpq)), aux$yzero)
  aux$zerosaux <- rep(0,max(aux$nmaxpq, aux$n+1))
  if(mean.correction){
    initial.values <- initial.values[-1]
    aux$upper <- aux$upper[-1]
    aux$lower <- aux$lower[-1]
  }
  aux$initial.values.arma <- initial.values
  if(is.null(objective.penalty)){
    aux$objective.penalty <- lgarchObjective(initial.values, aux)
  }else{
    aux$objective.penalty <- objective.penalty
  }

  #estimate:
  if(method=="ls"){
    objective.f <- function(pars, x=aux){ lgarchObjective(pars,x) }
  }else{
    objective.f <- function(pars, x=aux){ -lgarchObjective(pars,x) }
  }
  est <- nlminb(initial.values, objective.f, lower=aux$lower,
    upper=aux$upper, control=nlminb.control)

  #post-estimation:
  if(method!="ls") est$objective <- -est$objective
  if(mean.correction){
    est$par <- c(0,est$par)
    if(aux$ar > 0){
      arsum <- sum(est$par[aux$ar.indx])
    }else{
      arsum <- 0
    }
    est$par[1] <- (1-arsum)*aux$Elny2
  }
  names(est)[2] <- "objective.arma"

  #Elnz2 and sigma2u parameters:
  if(mean.correction){
    uadj <- lgarchRecursion1(c(0,est$par[-1]), aux)
  }else{
    uadj <- lgarchRecursion1(est$par, aux)
  }
  if(aux$yzeron > 0){
    uadj <- uadj[-aux$yzerowhere]
  }
  if(method=="ls"){
    Elnz2 <- -log(mean(exp(uadj)))
    namesArma <- NULL
    sigma2u <- var(uadj)
  }
  if(method=="ml"){
    Elnz2 <- -log(mean(exp(uadj - mean(uadj))))
    namesArma <- "sigma2u"
    sigma2u <- est$par[aux$sigma2u.indx]
  }
  if(method=="cex2"){
    Elnz2<- est$par[aux$sigma2u.indx]
    namesArma <- "Elnz2"
    sigma2u <- var(uadj)
  }
  par.lgarch <- Elnz2
  namesLgarch <- "Elnz2"
  par.arma <- est$par

  #xreg, garch, arch and intercept parameters:
  if(aux$xreg.k > 0){
    namesArma <- c(xreg.names, namesArma)
    namesLgarch <- c(xreg.names, namesLgarch)
    #OLD:
    #namesArma <- c(paste("xreg",1:aux$xreg.k,sep=""), namesArma)
    #namesLgarch <- c(paste("xreg",1:aux$xreg.k,sep=""), namesLgarch)
    par.lgarch <- c(est$par[aux$xreg.indx], par.lgarch)
  }
  if(aux$ma > 0){
    namesArma <- c(paste("ma",1:aux$ma,sep=""), namesArma)
    namesLgarch <- c(paste("garch",1:aux$ma,sep=""), namesLgarch)
    par.lgarch <- c(-est$par[aux$ma.indx], par.lgarch)
  }
  if(aux$ar > 0){
    namesArma <- c(paste("ar",1:aux$ar,sep=""), namesArma)
    namesLgarch <- c(paste("arch",1:aux$ar,sep=""), namesLgarch)
    par.tmp <- c(est$par[aux$ma.indx], rep(0, aux$ar-aux$ma))
    par.lgarch <- c(est$par[aux$ar.indx]+par.tmp, par.lgarch)
  }

  #intercept and completion:
  namesArma <- c("intercept", namesArma)
  #OLD: par.arma <- est$par
  names(par.arma) <- namesArma
  namesLgarch <- c("intercept", namesLgarch)
  const.lgarch <- est$par[1] - (1+sum(est$par[aux$ma.indx]))*Elnz2
  par.lgarch <- c(const.lgarch,par.lgarch)
  names(par.lgarch) <- namesLgarch
  est$par <- par.lgarch
  est <- c(list(date=date(), par.arma=par.arma), est)

  #vcov matrix:
  if(vcov){
    par.arma <- as.numeric(par.arma)
    if(mean.correction){
      par.arma <- par.arma[-1]
      namesArma <- namesArma[-1]
    }
    hessian.arma <- -optimHess(par.arma, objective.f)
    colnames(hessian.arma) <- namesArma
    rownames(hessian.arma) <- namesArma
    vcov.arma <- solve(-hessian.arma, tol=solve.tol)
    if(aux$method=="ls"){
      vcov.arma <- sigma2u*2*vcov.arma
    }
    est <- c(list(aux=aux, hessian.arma=hessian.arma,
      vcov.arma=vcov.arma, vcov.lgarch=NULL), est)
    est$vcov.lgarch <- vcov.lgarch(est, arma=FALSE)
  }else{
    est <- c(list(aux=aux, hessian.arma=NA,
      vcov.arma=NA, vcov.lgarch=NA), est)
  }

  #out:
  if(is.null(est$aux)){
    est <- c(list(aux=aux),est)
  }
  class(est) <- "lgarch"
  return(est)
}
