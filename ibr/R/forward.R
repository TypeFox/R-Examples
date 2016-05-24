forward <- function(formula,data,subset,criterion="gcv",df=1.5,Kmin=1,Kmax=1e+06,smoother="k",kernel="g",rank=NULL,control.par=list(),cv.options=list(),varcrit=criterion) {
  cl <- match.call() 
  mf <- match.call()
  m <- match(c("formula", "data", "subset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- "na.fail"
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  attr(mt,"intercept") <- 0
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) stop("offset are not allowed")
  if (is.empty.model(mt)) stop("no model ?")
  if (any(apply(mf,2,is.factor)))  stop("no factor are allowed")
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf)
  attributes(x) <- attributes(x)[c("dim","dimnames")]
  crit <-c("aic","aicc","gcv","bic","gmdl","rmse","map")
  varcrit <- match.arg(varcrit,crit)
  garde <- NULL
  ptot <- ncol(x)
  possible <- 1:ptot
  resultat <- matrix(Inf,ptot,ptot)
  colnames(resultat) <- colnames(x)
  letop <- Inf
  bandwidth <- NULL
  if (!((any(names(control.par)=="bandwidth")&&is.numeric(control.par$bandwidth))|((any(names(control.par)=="dftotal"))&&control.par$dftotal))&(smoother=="k"))
    bandwidth <- rep(0,ptot)
  for (iter in 1:ptot) {
    for (i in possible) {
      garden <- sort(c(garde,i))
      if ((iter>1)&&(!((any(names(control.par)=="bandwidth")&&is.numeric(control.par$bandwidth))|((any(names(control.par)=="dftotal"))&&control.par$dftotal))&((smoother=="k")&(kernel!="g")))) res3 <- ibr.fit(x[,garden,drop=FALSE],y,criterion,df,Kmin,Kmax,smoother,kernel,rank,control.par=c(control.par,list(bandwidth=bandwidth[garden])),cv.options)   else
      res3 <- ibr.fit(x[,garden,drop=FALSE],y,criterion,df,Kmin,Kmax,smoother,kernel,rank,control.par,cv.options)
      if ((iter==1)&&(!((any(names(control.par)=="bandwidth")&&is.numeric(control.par$bandwidth))|((any(names(control.par)=="dftotal"))&&control.par$dftotal))&((smoother=="k")&(kernel!="g"))))
        bandwidth[i] <- res3$bandwidth
      if (any(names(control.par)=="exhaustive")&&control.par$exhaustive) {
        allcrit <- switch(varcrit,rmse=res3$criteria$rmse,map=res3$criteria$map,aic=res3$criteria$aic,aicc=res3$criteria$aicc,gcv=res3$criteria$gcv,bic=res3$criteria$bic,gmdl=res3$criteria$gmdl)
        allcrit <- allcrit[res3$iter]
      } else {
        if (criterion==varcrit) allcrit <- res3$criteria else
        allcrit <- summary(res3,varcrit)$criteria
      }
      resultat[iter,i] <- allcrit
    }
    mini <- which.min(resultat[iter,])
    letopnew <-min(resultat[iter,])
    if (letopnew>letop) {
      resultat <- resultat[1:(iter-1),]
      if (iter>1) {
        class(resultat) <- c("forwardibr", "matrix")
      } else {
        class(resultat) <- c("forwardibr", "vector")
      }
      return(resultat) } else letop <- letopnew
    coordonnee <- which(possible==mini)
    garde <- sort(c(garde,possible[coordonnee]))
    possible <- possible[-coordonnee]
  }
  if (iter>1) {
    class(resultat) <- c("forwardibr", "matrix")
  } else {
    class(resultat) <- c("forwardibr", "vector")
  }
  return(resultat)
}

