heckitrob <-
function(outcome, selection, control=heckitrob.control())
{
  if (class(outcome) != "formula") {
    stop("argument 'outcome' must be a formula")
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
  m <- match(c("outcome", "data", "subset", "weights", "offset"), 
             names(mf), 0)
  mfO <- mf[c(1, m)]
  mfO$na.action <- na.pass
  mfO$drop.unused.levels <- TRUE
  mfO$na.action <- na.pass
  mfO[[1]] <- as.name("model.frame")
  names(mfO)[2] <- "formula"
  mfO <- eval(mfO, parent.frame())
  mtO <- attr(mfO, "terms")
  XO <- model.matrix(mtO, mfO)
  NXO <- ncol(XO)
  YO <- model.response(mfO)
  result$stage1 <- glmrob(YS~XS-1, family=binomial(link=probit),
                          method="Mqle",weghts.on.x=control$weights.x1, control = glmrobMqle.control(acc=control$acc, maxit=control$maxit, tcc=control$tcc))
  imrData=invMillsRatio(result$stage1)
  xMat=cbind(XO,imrData$IMR1)
  if(control$weights.x2=="none") x2weight=rep(1, length(YS)) else
    if(control$weights.x2=="hat") x2weight=sqrt(1-hat(xMat)) else
      if(control$weights.x2=="robCov") x2weight=x2weight.robCov(xMat) else
        if(control$weights.x2=="covMcd") x2weight=x2weight.covMcd(xMat)
  result$stage2=rlm(YO ~ XO+imrData$IMR1-1, method="M", psi=psi.huber, k=control$t.c, weights=x2weight, maxit=control$maxitO, subset=YS==1)
  xMat=model.matrix(result$stage2)
  x2weight=subset(x2weight, YS==1)
  result$vcov=heck2steprobVcov(YS[YS==1], YO[YS==1], model.matrix(result$stage1)[YS==1,], xMat, result$stage1, result$stage2$coeff, result$stage2$s, x2weight, control$t.c)
  result$method="robust two-stage"
  class(result)<- c("heckitrob", class(result))
  return(result)
}
