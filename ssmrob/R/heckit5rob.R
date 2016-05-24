heckit5rob <-
function(outcome1, outcome2, selection, control=heckitrob.control())
{
  if (class(outcome1) != "formula") {
    stop("argument 'outcome1' must be a formula")
  } 
  else  if (class(outcome2) != "formula") {
    stop("argument 'outcome2' must be a formula")
  }
  else if (length(selection) != 3) {
    stop("argument 'selection' must be a 2-sided formula")
  }  
  result=list()
    result$call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("selection", "data", "subset", "weights", "offset"), 
               names(mf), 0)
    mfS <- mf[c(1, m)]
    mfS$na.action <- na.pass
    mfS$drop.unused.levels <- TRUE
    mfS[[1]] <- as.name("model.frame")
    names(mfS)[2] <- "formula"
    mfS <- eval(mfS, parent.frame())
    mtS <- attr(mfS, "terms")
    XS <- model.matrix(mtS, mfS)
    NXS <- ncol(XS)
    YS <- model.response(mfS)
    m <- match(c("outcome1", "data", "subset", "weights", "offset"), 
               names(mf), 0)
    mfO1 <- mf[c(1, m)]
    mfO1$na.action <- na.pass
    mfO1$drop.unused.levels <- TRUE
    mfO1$na.action <- na.pass
    mfO1[[1]] <- as.name("model.frame")
    names(mfO1)[2] <- "formula"
    mfO1 <- eval(mfO1, parent.frame())
    mtO1 <- attr(mfO1, "terms")
    XO1 <- model.matrix(mtO1, mfO1)
    NXO1 <- ncol(XO1)
    YO1 <- model.response(mfO1)
  m <- match(c("outcome2", "data", "subset", "weights", "offset"), 
             names(mf), 0)
  mfO2 <- mf[c(1, m)]
  mfO2$na.action <- na.pass
  mfO2$drop.unused.levels <- TRUE
  mfO2$na.action <- na.pass
  mfO2[[1]] <- as.name("model.frame")
  names(mfO2)[2] <- "formula"
  mfO2 <- eval(mfO2, parent.frame())
  mtO2 <- attr(mfO2, "terms")
  XO2 <- model.matrix(mtO2, mfO2)
  NXO2 <- ncol(XO2)
  YO2 <- model.response(mfO2)
  result$stage1 <- glmrob(YS~XS-1, family=binomial(link=probit),
                            method="Mqle",weghts.on.x=control$weights.x1, control = glmrobMqle.control(acc=control$acc, maxit=control$maxit, tcc=control$tcc))
  imrData=invMillsRatio(result$stage1)
  IMR1=imrData$IMR1
  IMR2=-imrData$IMR0
  xMat1=cbind(XO1, IMR1)
  xMat2=cbind(XO2, IMR2)
  if(control$weights.x2=="none") x2weight1=rep(1, length(YS)) else
    if(control$weights.x2=="hat") x2weight1=sqrt(1-hat(xMat1)) else
      if(control$weights.x2=="robCov") x2weight1=x2weight.robCov(xMat1) else
        if(control$weights.x2=="covMcd") x2weight1=x2weight.covMcd(xMat1)
  if(control$weights.x2=="none") x2weight2=rep(1, length(YS)) else
    if(control$weights.x2=="hat") x2weight2=sqrt(1-hat(xMat2)) else
      if(control$weights.x2=="robCov") x2weight2=x2weight.robCov(xMat2) else
        if(control$weights.x2=="covMcd") x2weight2=x2weight.covMcd(xMat2)
  result$stage21=rlm(YO1 ~ XO1+IMR1-1, method="M", psi=psi.huber, k=control$t.c, weights=x2weight1, maxit=control$maxitO, subset=YS==1)
  result$stage22=rlm(YO2 ~ XO2+IMR2-1, method="M", psi=psi.huber, k=control$t.c, weights=x2weight2, maxit=control$maxitO, subset=YS==0)
  xMat1=model.matrix(result$stage21)
  xMat2=model.matrix(result$stage22)
  x2weight1=subset(x2weight1, YS==1)
  x2weight2=subset(x2weight2, YS==0)
  result$vcov1=heck2steprobVcov(YS[YS==1], YO1[YS==1], model.matrix(result$stage1)[YS==1,], xMat1, result$stage1, result$stage21$coeff, result$stage21$s, x2weight1, control$t.c)
  result$vcov2=heck5twosteprobVcov(YS[YS==0], YO2[YS==0], model.matrix(result$stage1)[YS==0,], xMat2, result$stage1, result$stage22$coeff, result$stage22$s, x2weight2, control$t.c)
  
  result$method="robust two-stage"
  class(result)<- c("heckit5rob", class(result))
  return(result)
}
