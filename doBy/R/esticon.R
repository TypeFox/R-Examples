esticon <- function(obj, cm, beta0, conf.int = TRUE, level=0.95, joint.test=FALSE,...)
  UseMethod("esticon")

esticon.gls <- function (obj, cm, beta0, conf.int = TRUE, level=0.95, joint.test=FALSE,...){
  if (joint.test==TRUE){
    .wald(obj, cm, beta0)
  } else {
    stat.name <- "X2.stat"
    vcv <- vcov(obj)
    coef.mat  <- matrix(coef(obj))
    df  <- 1
    .esticonCore(obj, cm, beta0, conf.int=conf.int,level,coef.mat,vcv,df,stat.name)
  }
}

esticon.geeglm <- function (obj, cm, beta0, conf.int = TRUE, level=0.95, joint.test=FALSE,...){
  if (joint.test==TRUE){
    .wald(obj, cm, beta0)
  } else {
    stat.name <- "X2.stat"
    coef.mat  <- summary(obj)$coef
    vcv <- summary(obj)$cov.scaled
    df  <- 1
    .esticonCore(obj, cm, beta0, conf.int=conf.int,level,coef.mat,vcv,df,stat.name)
  }
}

esticon.lm <- function (obj, cm, beta0, conf.int = TRUE, level=0.95, joint.test=FALSE,...){
  if (joint.test==TRUE){
    .wald(obj, cm, beta0)
  } else {
    stat.name <- "t.stat"
    coef.mat  <- summary.lm(obj)$coefficients
    coef.vec  <- coef(obj)
    vcv <- summary.lm(obj)$cov.unscaled * summary.lm(obj)$sigma^2
    df  <- obj$df.residual
    .esticonCore(obj, cm, beta0, conf.int=conf.int, level, coef.mat, vcv, df, stat.name, coef.vec=coef.vec)
  }
}

esticon.glm <- function (obj, cm, beta0, conf.int = TRUE, level=0.95, joint.test=FALSE,...){
  if (joint.test==TRUE){
    .wald(obj, cm, beta0)
  } else {
    coef.mat <- summary.lm(obj)$coefficients
    vcv <- summary(obj)$cov.scaled
    if(family(obj)[1] %in% c("poisson","binomial")){
      stat.name <- "X2.stat"
      df <- 1
    } else {
      stat.name <- "t.stat"
      df <- obj$df.residual
    }
    .esticonCore(obj, cm, beta0, conf.int=conf.int,level=level,coef.mat,vcv,df,stat.name)
  }
}

esticon.mer <- esticon.merMod <- function (obj, cm, beta0, conf.int = TRUE, level=0.95, joint.test=FALSE,...){
  if (joint.test==TRUE){
    .wald(obj, cm, beta0)
  } else {
    stat.name <- "X2.stat"
    coef.mat  <- matrix(lme4::fixef(obj))
    vcv <- as.matrix(vcov(obj))
    df  <- 1
    .esticonCore(obj, cm, beta0, conf.int=conf.int, level, coef.mat, vcv, df, stat.name)
  }
}

esticon.coxph <-
  function (obj, cm, beta0, conf.int = TRUE, level = 0.95, joint.test = FALSE, ...)
{
  if (joint.test == TRUE) {
    .wald(obj, cm, beta0)
  }
  else {
    cf <- summary(obj)$coefficients
    vcv <- obj$var
    stat.name <- "X2.stat"
    df <- 1
    .esticonCore(obj, cm, beta0, conf.int = conf.int, level = level,
                 cf, vcv, df, stat.name)
  }
}



### ######################################################
###
### esticon.lme needs better testing at some point
###
### ######################################################

esticon.lme <- function (obj, cm, beta0, conf.int = NULL, level=0.95, joint.test=FALSE,...){
  warning("The esticon function has not been thoroughly teste on 'lme' objects")
  if (joint.test==TRUE){
    .wald(obj, cm, beta0)
  } else {
    stat.name <- "t.stat"
    coef.mat <- summary(obj)$tTable
    rho <- summary(obj)$cor
    vcv <- rho * outer(coef.mat[, 2], coef.mat[, 2])
    tmp <- cm
    tmp[tmp == 0] <- NA
    df.all <- t(abs(t(tmp) * obj$fixDF$X))
    df <- apply(df.all, 1, min, na.rm = TRUE)
    problem <- apply(df.all != df, 1, any, na.rm = TRUE)
    if (any(problem))
      warning(paste("Degrees of freedom vary among parameters used to ",
                    "construct linear contrast(s): ",
                    paste((1:nrow(tmp))[problem],
                          collapse = ","),
                    ". Using the smallest df among the set of parameters.",
                    sep = ""))
    df <- min(df)
    .esticonCore(obj, cm, beta0, conf.int=conf.int,level,coef.mat,vcv,df,stat.name)
  }
}


### .functions below here ###

.wald <- function (obj, cm,beta0)
{
    if (!is.matrix(cm) && !is.data.frame(cm))
        cm <- matrix(cm, nrow = 1)

    if (missing(beta0))
      beta0 <- rep(0,nrow(cm))

    df <- nrow(cm)
    if ("geese" %in% class(obj)) {
      coef.mat  <- obj$beta
      vcv <- obj$vbeta
    } else if ("geeglm" %in% class(obj)) {
      coef.mat  <- obj$coef
      vcv <- summary(obj)$cov.scaled
    } else if ("gls" %in% class(obj)) {
      vcv <- vcov(obj)
      coef.mat  <- matrix(coef(obj))
    } else if ("gee" %in% class(obj)) {
      coef.mat  <- obj$coef
      vcv <- obj$robust.variance
    }
    else if ("lm" %in% class(obj)) {
      coef.mat  <- summary.lm(obj)$coefficients[, 1]
      vcv <- summary.lm(obj)$cov.unscaled * summary.lm(obj)$sigma^2
      if ("glm" %in% class(obj)) {
        vcv <- summary(obj)$cov.scaled
      }
    }
    else if ("coxph" %in% class(obj)) {
      coef.mat <- obj$coef
      vcv <- obj$var
    }
    else
      stop("obj must be of class 'lm', 'glm', 'aov', 'gls', 'gee', 'geese', 'coxph'")
    u      <- (cm %*% coef.mat)-beta0
    vcv.u  <- cm %*% vcv %*% t(cm)
    W      <- t(u) %*% solve(vcv.u) %*% u
    prob   <- 1 - pchisq(W, df = df)
    retval <- as.data.frame(cbind(W, df, prob))
    names(retval) <- c("X2.stat", "DF", "Pr(>|X^2|)")
    return(as.data.frame(retval))
}


.esticonCore <- function (obj, cm, beta0, conf.int = NULL, level,coef.mat ,vcv,df,stat.name, coef.vec=coef.mat[,1] ) {

  ## Notice
  ## coef.mat: summary(obj)$coefficients
  ## coef.vec: coef(obj)
  ## cl <- match.call(); print(cl);  print(coef.mat)

  if (missing(cm)) stop("Contrast matrix 'cm' is missing")
  if (!is.matrix(cm) && !is.data.frame(cm))
    cm <- matrix(cm, nrow = 1)
  if (missing(beta0))
    beta0 <- rep(0,nrow(cm))

  idx <- !is.na(coef.vec) ## Only want columns of cm for identifiable parameters
  cm <- cm[ , idx, drop=FALSE]

  if (!dim(cm)[2] == dim(coef.mat)[1])
    stop(paste("\n Dimension of ",
               deparse(substitute(cm)),
               ": ", paste(dim(cm), collapse = "x"),
               ", not compatible with no of parameters in ",
               deparse(substitute(obj)), ": ", dim(coef.mat)[1], sep = ""))
  ct      <- cm %*% coef.mat[, 1]
  ct.diff <- cm %*% coef.mat[, 1] - beta0
  vc      <- sqrt(diag(cm %*% vcv %*% t(cm)))

  switch(stat.name,
         t.stat = {
           prob <- 2 * (1 - pt(abs(ct.diff/vc), df))
         },
         X2.stat = {
           prob <- 1 - pchisq((ct.diff/vc)^2, df = 1)
         })

  if (stat.name == "X2.stat") {
    retval <- cbind(hyp=beta0, est = ct, stderr = vc, t.value = (ct.diff/vc)^2,
                    df = df, prob = prob )
    dimnames(retval) <-
      list(NULL, c("beta0","Estimate","Std.Error","X2.value","DF","Pr(>|X^2|)"))
  }
  else if (stat.name == "t.stat") {
    retval <- cbind(hyp=beta0, est = ct, stderr = vc, t.value = ct.diff/vc,
                    df = df, prob = prob)
    dimnames(retval) <-
      list(NULL, c("beta0","Estimate","Std.Error","t.value","DF","Pr(>|t|)"))
  }

  conf.int <- "wald"
  if (!is.null(conf.int)) {
    if (level <= 0 || level >= 1)
      stop("level should be betweeon 0 and 1. Usual values are 0.95, 0.90")

    alpha <- 1 - level
    switch(stat.name,
           t.stat  = { quant <- qt(1 - alpha/2, df  )  },
           X2.stat = { quant <- qnorm(1 - alpha/2) })

    vvv <- cbind(ct.diff-vc*quant, ct.diff+vc*quant)
    colnames(vvv) <- c("Lower", "Upper")
    retval <- cbind(retval, vvv)
  }
  retval[,6] <- round(retval[,6],7)
  return(as.data.frame(retval))
}


##  rn <- NULL
##   if (is.null(rownames(cm)))
##     rn <- paste("(", apply(cm, 1, paste, collapse = " "), ")", sep = "")
##   else
##     rn <- rownames(cm)








