predict.spatgev <- function(object, newdata, ret.per = NULL,
                            ...){

  param <- object$param
    
  if (!missing(newdata)){
    data <- newdata
    marg.cov <- NULL
    
    if (is.null(dim(newdata)))
      data <- t(as.matrix(data))
  }

  else
    data <- object$covariables

  loc.form <- object$loc.form
  scale.form <- object$scale.form
  shape.form <- object$shape.form
  
  loc.dsgnmat <- modeldef(data, loc.form)$dsgn.mat
  scale.dsgnmat <- modeldef(data, scale.form)$dsgn.mat
  shape.dsgnmat <- modeldef(data, shape.form)$dsgn.mat
  
  idx.loc <- which(substr(names(param), 1, 8) == "locCoeff")
  idx.scale <- which(substr(names(param), 1, 10) == "scaleCoeff")
  idx.shape <- which(substr(names(param), 1, 10) == "shapeCoeff")
  
  loc.pred <- loc.dsgnmat %*% param[idx.loc]
  scale.pred <- scale.dsgnmat %*% param[idx.scale]
  shape.pred <- shape.dsgnmat %*% param[idx.shape]

  ans <- cbind(loc.pred, scale.pred, shape.pred)
  colnames(ans) <- c("loc", "scale", "shape")
  ans <- cbind(data, ans)

  if (!is.null(ret.per)){
    ret.lev <- NULL
    for (T in ret.per)
      ret.lev <- cbind(ret.lev, .qgev(1 - 1/T, ans[,"loc"],
                                      ans[,"scale"],
                                      ans[,"shape"]))

  colnames(ret.lev) <- paste("Q", ret.per, sep="")
  ans <- cbind(ans, ret.lev)
  }
  
  return(ans)
}

predict.maxstab <- function(object, newdata, ret.per = NULL,
                            std.err = TRUE, ...){
  
  param <- object$param
    
  if (!missing(newdata)){
    data <- newdata
    marg.cov <- NULL
    
    if (is.null(dim(newdata)))
      data <- t(as.matrix(data))
  }

  else{
    data <- object$coord
    marg.cov <- object$marg.cov
  }

  if (!is.null(marg.cov))
    data <- cbind(data, marg.cov)

  if (object$fit.marge){
    loc.form <- object$loc.form
    scale.form <- object$scale.form
    shape.form <- object$shape.form
    
    loc.dsgnmat <- modeldef(data, loc.form)$dsgn.mat
    scale.dsgnmat <- modeldef(data, scale.form)$dsgn.mat
    shape.dsgnmat <- modeldef(data, shape.form)$dsgn.mat
    
    idx.loc <- which(substr(names(param), 1, 8) == "locCoeff")
    idx.scale <- which(substr(names(param), 1, 10) == "scaleCoeff")
    idx.shape <- which(substr(names(param), 1, 10) == "shapeCoeff")
    
    loc.pred <- loc.dsgnmat %*% param[idx.loc]
    scale.pred <- scale.dsgnmat %*% param[idx.scale]
    shape.pred <- shape.dsgnmat %*% param[idx.shape]

    if (std.err && is.matrix(object$var.cov)){
      par.est <- object$fitted
      idx.loc <- which(substr(names(par.est), 1, 8) == "locCoeff")
      idx.scale <- which(substr(names(par.est), 1, 10) == "scaleCoeff")
      idx.shape <- which(substr(names(par.est), 1, 10) == "shapeCoeff")

      loc.std.err <- sqrt(diag(loc.dsgnmat %*% object$var.cov[idx.loc, idx.loc] %*%
                               t(loc.dsgnmat)))
      scale.std.err <- sqrt(diag(scale.dsgnmat %*% object$var.cov[idx.scale, idx.scale] %*%
                                 t(scale.dsgnmat)))
      shape.std.err <- sqrt(diag(shape.dsgnmat %*% object$var.cov[idx.shape, idx.shape] %*%
                                 t(shape.dsgnmat)))
    }

    else
      loc.std.err <- scale.std.err <- shape.std.err <- NULL
  }

  else {
    loc.pred <- scale.pred <- shape.pred <- rep(1, nrow(data))
    loc.std.err <- scale.std.err <- shape.std.err <- NULL
  }

  ans <- cbind(loc.pred, scale.pred, shape.pred)
  colnames(ans) <- c("loc", "scale", "shape")
  ans <- cbind(data, ans)

  if (!is.null(loc.std.err))
    ans <- cbind(ans, loc.std.err = loc.std.err, scale.std.err = scale.std.err,
                 shape.std.err = shape.std.err)

  
  if (!is.null(ret.per)){
    ret.lev <- NULL
    for (T in ret.per)
      ret.lev <- cbind(ret.lev, .qgev(1 - 1/T, ans[,"loc"],
                                      ans[,"scale"],
                                      ans[,"shape"]))

    colnames(ret.lev) <- paste("Q", ret.per, sep="")
    ans <- cbind(ans, ret.lev)

    if (std.err && is.matrix(object$var.cov)){
      par.est <- object$fitted
      idx.loc <- which(substr(names(par.est), 1, 8) == "locCoeff")
      idx.scale <- which(substr(names(par.est), 1, 10) == "scaleCoeff")
      idx.shape <- which(substr(names(par.est), 1, 10) == "shapeCoeff")
      idx.all <- c(idx.loc, idx.scale, idx.shape)

      ret.lev.std.err <- NULL
      for (i in 1:length(ret.per)){
        gradient <- cbind(loc.dsgnmat, (ret.lev[,i] - ans[,"loc"]) / ans[,"scale"] *
                          scale.dsgnmat,
                          (-log(-log(1 - 1 /ret.per[i])) * ans[,"scale"] / ans[,"shape"] *
                         (- log(1 - 1 / ret.per[i]))^(-ans[,"shape"]) -
                         (ret.lev[,i] - ans[,"loc"]) / ans[,"shape"]) * shape.dsgnmat)

      ret.lev.std.err <- cbind(ret.lev.std.err,
                               sqrt(diag(gradient %*% object$var.cov[idx.all, idx.all] %*%
                                         t(gradient))))
      }

      colnames(ret.lev.std.err) <- paste("Q", ret.per, "std.err", sep ="")
      ans <- cbind(ans, ret.lev.std.err)
    }
    
    else
      ret.lev.std.err <- NULL
    }
  
  return(ans)
}


predict.copula <- function(object, newdata, ret.per = NULL,
                           std.err = TRUE, ...){
  
  param <- object$param
  
  if (!missing(newdata)){
    data <- newdata
    marg.cov <- NULL
    
    if (is.null(dim(newdata)))
      data <- t(as.matrix(data))
  }
  
  else{
    data <- object$coord
    marg.cov <- object$marg.cov
  }
  
  if (!is.null(marg.cov))
    data <- cbind(data, marg.cov)
  
  if (object$fit.marge){
    loc.form <- object$loc.form
    scale.form <- object$scale.form
    shape.form <- object$shape.form
    
    loc.dsgnmat <- modeldef(data, loc.form)$dsgn.mat
    scale.dsgnmat <- modeldef(data, scale.form)$dsgn.mat
    shape.dsgnmat <- modeldef(data, shape.form)$dsgn.mat
    
    idx.loc <- which(substr(names(param), 1, 8) == "locCoeff")
    idx.scale <- which(substr(names(param), 1, 10) == "scaleCoeff")
    idx.shape <- which(substr(names(param), 1, 10) == "shapeCoeff")
    
    loc.pred <- loc.dsgnmat %*% param[idx.loc]
    scale.pred <- scale.dsgnmat %*% param[idx.scale]
    shape.pred <- shape.dsgnmat %*% param[idx.shape]
    
    if (std.err && is.matrix(object$var.cov)){
      par.est <- object$fitted
      idx.loc <- which(substr(names(par.est), 1, 8) == "locCoeff")
      idx.scale <- which(substr(names(par.est), 1, 10) == "scaleCoeff")
      idx.shape <- which(substr(names(par.est), 1, 10) == "shapeCoeff")
      
      loc.std.err <- sqrt(diag(loc.dsgnmat %*% object$var.cov[idx.loc, idx.loc] %*%
                               t(loc.dsgnmat)))
      scale.std.err <- sqrt(diag(scale.dsgnmat %*% object$var.cov[idx.scale, idx.scale] %*%
                                 t(scale.dsgnmat)))
      shape.std.err <- sqrt(diag(shape.dsgnmat %*% object$var.cov[idx.shape, idx.shape] %*%
                                 t(shape.dsgnmat)))
    }
    
    else
      loc.std.err <- scale.std.err <- shape.std.err <- NULL
  }
  
  else {
    loc.pred <- scale.pred <- shape.pred <- rep(1, nrow(data))
    loc.std.err <- scale.std.err <- shape.std.err <- NULL
  }

  ans <- cbind(loc.pred, scale.pred, shape.pred)
  colnames(ans) <- c("loc", "scale", "shape")
  ans <- cbind(data, ans)

  if (!is.null(loc.std.err))
    ans <- cbind(ans, loc.std.err = loc.std.err, scale.std.err = scale.std.err,
                 shape.std.err = shape.std.err)
  
  if (!is.null(ret.per)){
    ret.lev <- NULL
    for (T in ret.per)
      ret.lev <- cbind(ret.lev, .qgev(1 - 1/T, ans[,"loc"],
                                      ans[,"scale"],
                                      ans[,"shape"]))
    
    colnames(ret.lev) <- paste("Q", ret.per, sep="")
    ans <- cbind(ans, ret.lev)
  }
  
  return(ans)
}

predict.pspline <- function(object, newdata, ...){

  if (missing(newdata))
    new.data <- object$x

  degree <- object$degree
  knots <- object$knots
  beta <- object$beta
  
  dsgn.mat <- rb(newdata, degree = degree, knots = knots,
                 penalty = NULL)$dsgn.mat

  y <- dsgn.mat %*% beta
  ans <- cbind(newdata, y)
  colnames(ans) <- c("x", "y.hat")
  
  return(ans)
}
