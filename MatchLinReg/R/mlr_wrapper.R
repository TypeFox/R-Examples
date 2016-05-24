colname.dissemble <- function(str) strsplit(str, split = ":")[[1]]
colname.assemble <- function(str.vec) paste(str.vec, collapse = ":")
colname.standardize <- function(str) colname.assemble(sort(colname.dissemble(str)))

mlr.match <- function(tr, X, psm = TRUE, replace = F, caliper = Inf, verbose = TRUE) { # consider removing ... to tighten control (e.g. no replacement)
  X <- as.matrix(X)
  if (psm) {
    psm.reg <- glm(tr ~ X, family = "binomial")
    my.match <- Match(Tr = tr, X = psm.reg$linear.predictors, replace = replace, caliper = caliper, M = 1)
  } else {
    psm.reg <- NA
    my.match <- Match(Tr = tr, X = X, replace = replace, caliper = caliper, M = 1)
  }
  ret <- c(my.match$index.control, my.match$index.treated)
  
  nc <- length(my.match$index.control)
  nt <- length(my.match$index.treated)
  
  if (verbose) {
    cat("size of control group:", nc, "(", 100 * nc / length(which(tr ==0)), " %)\n")
    cat("size of treatment group:", nt, "(", 100 * nt / length(which(tr == 1)), " %)\n")
  }
  
  # add details flag
  attr(ret, "nc") <- nc
  attr(ret, "nt") <- nt
  attr(ret, "psm.reg") <- psm.reg
  attr(ret, "match.obj") <- my.match
  
  return (ret)
}

# 1: relative squared bias reduction by term, at best index (for what omitted R-squared value?)
# 2: squared bias by term vs. index
# 3: combined bias (all three methods) vs. index
# 4: bias, variance and MSE vs. index, for multiple values of omitted R-squared
# 5: best index vs. omitted R -squared
summary.mlr <- function(object, power = FALSE
  , power.control = list(rnd = TRUE, d = 0.5, sig.level = 0.05, niter = 1000, rnd = TRUE)
  , max.method = c("single-covariate", "covariate-subspace", "absolute")
  , verbose = FALSE, ...
  , orsq.min = 1e-03, orsq.max = 1e0, n.orsq = 100) {
  tr <- object$tr
  Z.i <- object$Z.i
  Z.o <- object$Z.o
  idx.list <- object$idx.list
  nidx <- length(idx.list)
  nterms <- ncol(Z.o)
  
  Z.o.orth <- mlr.orthogonalize(X = Z.i, Z = Z.o, normalize = T)
  bias.obj.baseline <- mlr.bias(tr = tr, Z.i = Z.i, Z.o = Z.o.orth)
  bias.baseline <- c(bias.obj.baseline$single$bias, bias.obj.baseline$subspace$bias
                     , bias.obj.baseline$absolute$bias)
  bias.terms.baseline <- bias.obj.baseline$single$bias.vec
  Z.single <- bias.obj.baseline$single$dir
  Z.subspace <- bias.obj.baseline$subspace$dir
  Z.absolute <- bias.obj.baseline$absolute$dir
  var.baseline <- mlr.variance(tr = tr, Z.i = Z.i, sigsq = 1, details = F)
  smd.baseline <- mlr.smd(tr, cbind(Z.i, Z.o))
  
  bias.mat <- array(NA, dim = c(nidx, 3))
  bias.terms.mat <- array(NA, dim = c(nidx, nterms))
  var.vec <- rep(NA, nidx)
  nt.vec <- rep(NA, nidx)
  nc.vec <- rep(NA, nidx)
  smd.mat <- array(NA, dim = c(nidx, ncol(Z.i) + ncol(Z.o)))
  rettmp <- sapply(1:nidx, function(n) {
    idx <- idx.list[[n]]
    bias.mat[n, ] <<- c(mlr.bias(tr = tr[idx], Z.i = Z.i[idx, ], Z.o = Z.single[idx], gamma.o = 1.0)$gamma.o,
                           mlr.bias(tr = tr[idx], Z.i = Z.i[idx, ], Z.o = Z.subspace[idx], gamma.o = 1.0)$gamma.o,
                           mlr.bias(tr = tr[idx], Z.i = Z.i[idx, ], Z.o = Z.absolute[idx], gamma.o = 1.0)$gamma.o)
    bias.terms.mat[n, ] <<- mlr.bias(tr = tr[idx], Z.i = Z.i[idx, ], Z.o = Z.o.orth[idx, ], gamma.o = rep(1, nterms))$single$bias.vec
    var.vec[n] <<- mlr.variance(tr = tr[idx], Z.i = Z.i[idx, ], sigsq = 1.0, details = F)
    nt.vec[n] <<- length(which(tr[idx] == 1))
    nc.vec[n] <<- length(which(tr[idx] == 0))
    smd.mat[n, ] <<- mlr.smd(tr[idx], cbind(Z.i, Z.o)[idx, ])
    
    if (verbose) {
      cat("treatment/control group sizes: ", nt.vec[n], "/", nc.vec[n], "\n", sep = "")
    }
    
    return (NULL)
  })
  bias.mat <- rbind(bias.mat, bias.baseline)
  rownames(bias.mat) <- NULL; colnames(bias.mat) <- c("single-covariate", "covariate-subspace", "absolute")
  bias.terms.mat <- rbind(bias.terms.mat, bias.terms.baseline)
  rownames(bias.terms.mat) <- NULL
  var.vec <- c(var.vec, var.baseline)
  smd.mat <- rbind(smd.mat, smd.baseline)
  rownames(smd.mat) <- NULL
  
  if (power) {
    power.baseline <- mlr.power(tr = tr, Z.i = Z.i, d = power.control$d, sig.level = power.control$sig.level, niter = power.control$niter, rnd = power.control$rnd)
    power.mat <- array(NA, dim = c(nidx, length(power.baseline)))
    rettmp <- sapply(1:nidx, function(n) {
      idx <- idx.list[[n]]
      power.mat[n, ] <<- mlr.power(tr = tr, Z.i = Z.i, d = power.control$d, sig.level = power.control$sig.level, niter = power.control$niter, rnd = power.control$rnd, idx = idx)
    })
    power.mat <- rbind(power.mat, power.baseline)
    colnames(power.mat) <- c("mean", "se", "rnd.mean", "rnd.se")
    rownames(power.mat) <- NULL
  } else {
    power.mat <- NA
  }
  
  # combining bias and variance for each value of omitted r-squared
  max.method <- match.arg(max.method)
  if (max.method == "single-covariate") {
    max.index <- 1
  } else if (max.method == "covariate-subspace") {
    max.index <- 2
  } else if (max.method == "absolute") {
    max.index <- 3
  }
  bvmat <- cbind(bias.mat[, max.index], var.vec)
  combine.obj <- mlr.combine.bias.variance(tr = tr, bvmat = bvmat, orsq.min = orsq.min, orsq.max = orsq.max, n.orsq = n.orsq) # create parameters
  
  ret <- list(mlr.obj = object, bias = bias.mat, bias.terms = bias.terms.mat, variance = var.vec, power = power.mat, smd = smd.mat
              , combine.obj = combine.obj)
  class(ret) <- "summary.mlr"
  
  return (ret)
}

plot.summary.mlr <- function(x, which = 1
                             , smd.index = 1:min(10, ncol(x$smd))
                             , bias.index = 1:min(10, ncol(x$bias.terms))
                             , orsq.plot = c(0.01, 0.05, 0.25)
                             , caption.vec = c("relative squared bias reduction", "normalized bias"
                                               , "standardized mean difference", "maximum bias"
                                               , "error components", "optimum choice", "power analysis")
                             , ...) {
  #delta.x.1 <- 3.0
  #delta.x.2 <- 1 * delta.x.1
  #xvec.orig <- x$mlr.obj$caliper.vec
  #xvec <- xvec.orig
  #xvec[is.infinite(xvec)] <- max(xvec[!is.infinite(xvec)]) * delta.x.1 # change this so xvec is produced by mlr not here
  #xvec <- c(xvec, max(xvec) * delta.x.2) # for no matching
  #xlab <- format(xvec.orig, digits = 2); xlab[is.infinite(xvec.orig)] <- "no caliper"; xlab <- c(xlab, "no match")
  #xtitle <- "caliper size" # this should also come from mlr
  
  xvec <- x$mlr.obj$xvec
  xlab <- x$mlr.obj$xlab
  xtitle <- x$mlr.obj$xtitle
  
  cex.mult <- 1.3
  caption.vec <- rep(caption.vec, 7)
  
  which.list <- 1:7
  which.valid <- which.list[which.list %in% which]
  which.rsbr <- 1 # relative squared bias reduction (by term) for a single idx
  which.terms <- 2 # bias terms vs. idx
  which.smd <- 3 # standardized mean difference by term vs. idx
  which.max <- 4 # maximum bias (single-covariate, subspace, absolute) vs. idx
  which.combine <- 5 # bias/variance/MSE plots
  which.calibrate <- 6 # optimum index vs. omitted r-squared
  which.power <- 7 # power analysis (matched and random subsamples) vs. idx
  
  # optimum index vs. omitted r-squared
  if (which.calibrate %in% which.valid) {
    combine.obj <- x$combine.obj
    plot(combine.obj$orsq.vec, xvec[combine.obj$which.min.vec], type = "l", log = "xy", xlab = "omitted r-squared"
         , ylab = paste("optimal ", xtitle, sep = ""), yaxt = "n"
         , main = caption.vec[which.calibrate]
         , cex = cex.mult, cex.axis = cex.mult, cex.lab = cex.mult, cex.main = cex.mult
         )
    axis(side = 2, at = xvec, labels = xlab
         #, cex.axis = cex.mult
         )
  }
  
  # bias/variance/MSE plots
  if (which.combine %in% which.valid) {
    pch.vec <- 1:3
    combine.obj <- x$combine.obj
    orsq.index.vec <- sapply(orsq.plot, function(orsq) {
      abs.dist <- abs(combine.obj$orsq.vec - orsq)
      return (which(abs.dist == min(abs.dist))[1])
    })
    
    for (orsq.index in orsq.index.vec) {
      plot(xvec, combine.obj$errmat[orsq.index, ], type = "b"
           , ylim = range(combine.obj$errmat[orsq.index, ], x$variance, combine.obj$biassq.mat[orsq.index, ])
           , xlab = xtitle, ylab = "TE error", xaxt = "n"
           , main = paste(caption.vec[which.combine], ", o-rsq = ", format(combine.obj$orsq.vec[orsq.index], digits = 2), sep = "")
           , pch = pch.vec[1]
           , cex = cex.mult, cex.axis = cex.mult, cex.lab = cex.mult, cex.main = cex.mult
           , log = "x"
           )
      axis(side = 1, at = xvec, labels = xlab)
      lines(xvec, x$variance, type = "b", pch = pch.vec[2])
      lines(xvec, combine.obj$biassq.mat[orsq.index, ], type = "b", pch = pch.vec[3])
      legend("left", legend = c("total MSE", "variance", "bias"), pch = pch.vec, lty = rep(1, 3), cex = cex.mult)
    }
  }
  
  # maximum bias (single-covariate, subspace, absolute) vs. idx
  if (which.max %in% which.valid) {
    bias.mat.sq <- x$bias^2
    pch.vec <- 1:3
    plot(xvec, bias.mat.sq[, 1], type = "b", xlab = xtitle, ylab = "normalized squared bias", xaxt = "n", log = "xy"
         , pch = pch.vec[1], ylim = range(bias.mat.sq)
         , cex = cex.mult, cex.axis = cex.mult, cex.lab = cex.mult, cex.main = cex.mult
         , main = caption.vec[which.max]
         )
    axis(side = 1, at = xvec, labels = xlab)
    lines(xvec, bias.mat.sq[, 2], type = "b", pch = pch.vec[2])
    lines(xvec, bias.mat.sq[, 3], type = "b", pch = pch.vec[3])
    legend("bottomright", legend = colnames(bias.mat.sq), pch = pch.vec, lty = rep(1,3))
  }
  
  # standardized mean difference, by terms vs. idx
  if (which.smd %in% which.valid) {
    smd.index <- smd.index[smd.index %in% 1:ncol(x$smd)]
    smd.mat.abs <- abs(x$smd[, smd.index, drop = F])
    pch.vec <- 1:ncol(smd.mat.abs)
    lty.vec <- rep(1, ncol(smd.mat.abs))
    plot(xvec, smd.mat.abs[, 1], type = "b", xlab = xtitle, xaxt = "n", ylab = "absolute smd"
         #, xlim = c(0.0, 2.0)
         , ylim = range(smd.mat.abs), pch = pch.vec[1], log = "x"
         , cex = cex.mult, cex.axis = cex.mult, cex.lab = cex.mult, cex.main = cex.mult
         , main = caption.vec[which.smd]
         )
    axis(side = 1, at = xvec, labels = xlab)
    if (ncol(smd.mat.abs) > 1) {
      for (i in 2:ncol(smd.mat.abs)) {
        lines(xvec, smd.mat.abs[, i], type = "b", pch = pch.vec[i])
      }
    }
    legend("topleft", legend = colnames(smd.mat.abs), lty = lty.vec, pch = pch.vec)
  }
  
  # relative squared bias reduction
  if (which.rsbr %in% which.valid) {
    bias.terms <- x$bias.terms
    rsbr.idx <- nrow(bias.terms) - 1 # obtain this after optimization, or from user
    bias.terms.nomatch <- bias.terms[nrow(bias.terms), ]
    bias.terms.match <- bias.terms[rsbr.idx, ]
    rsbr <- (bias.terms.nomatch^2 - bias.terms.match^2) / bias.terms.nomatch^2
    rsbr <- sort(rsbr, decreasing = T)
    which.rsbr.improve <- which(rsbr > 0)
    which.rsbr.worse <- which(rsbr <= 0)
    improve.range <- range(rsbr[which.rsbr.improve])
    yrange <- c(improve.range[1] - 0.4, +1)
    text.y <- improve.range[1] - 0.2
    
    plot(rsbr[which.rsbr.improve], xaxt = "n", xlab = "candidate omitted term", ylab = "relative squared bias reduction"
         , ylim = yrange, type = "b", pch = 3
         , xlim = c(1, length(rsbr))
         , cex = cex.mult, cex.axis = cex.mult, cex.lab = cex.mult, cex.main = cex.mult
         , main = caption.vec[which.rsbr]
         )
    #axis(1, at = 1:length(which.improve), labels = FALSE)
    text(x = which.rsbr.improve, y = text.y, labels = names(rsbr)[which.rsbr.improve]
         , srt = 90, pos = 1, xpd = T)
    text(x = which.rsbr.worse, y = text.y, labels = names(rsbr)[which.rsbr.worse]
         , srt = 90, pos = 1, xpd = T, col = "red")
    text(x = which.rsbr.worse, y = min(rsbr[which.rsbr.improve])
         , labels = format(rsbr[which.rsbr.worse], digits = 1), srt = 90, pos = 1, xpd = T, col = "red")
  }
  
  # bias terms
  if (which.terms %in% which.valid) {
    bias.index <- bias.index[bias.index %in% 1:ncol(x$bias.terms)]
    bias.sq.terms <- (x$bias.terms[, bias.index, drop = F])^2
    pch.vec <- 1:ncol(bias.sq.terms)
    lty.vec <- rep(1, ncol(bias.sq.terms))
    plot(xvec, bias.sq.terms[, 1], type = "b", xlab = xtitle, xaxt = "n", ylab = "normalized squared bias"
         #, xlim = c(0.0, 2.0)
         , ylim = range(bias.sq.terms), pch = pch.vec[1], log = "x"
         , cex = cex.mult, cex.axis = cex.mult, cex.lab = cex.mult, cex.main = cex.mult
         , main = caption.vec[which.terms]
         )
    axis(side = 1, at = xvec, labels = xlab)
    if (ncol(bias.sq.terms) > 1) {
      for (i in 2:ncol(bias.sq.terms)) {
        lines(xvec, bias.sq.terms[, i], type = "b", pch = pch.vec[i])
      }
    }
    legend("topleft", legend = colnames(bias.sq.terms), lty = lty.vec, pch = pch.vec)
  }

  # power analysis
  if (which.power %in% which.valid && !is.na(x$power)) {
    power.mat <- x$power
    power.mean <- power.mat[, 1]
    power.rnd.mean <- power.mat[, 3]
    power.se <- power.mat[, 2]
    power.rnd.se <- power.mat[, 4]
    power.plus <- power.mean + 2 * power.se
    power.minus <- power.mean - 2 * power.se
    power.rnd.plus <- power.rnd.mean + 2 * power.rnd.se
    power.rnd.minus <- power.rnd.mean - 2 * power.rnd.se
    errbar(x = xvec, y = power.mean, yplus = power.plus, yminus = power.minus
           , xaxt = "n", xlab = xtitle, ylab = "power", log = "x"
           , cex = cex.mult, cex.axis = cex.mult, cex.lab = cex.mult, cex.main = cex.mult
           , lty = 1, type = "b"
           , ylim = range(power.mat[, c(1,3)])
           #, main = caption.vec[which.power]
           )
    title(main = caption.vec[which.power])
    axis(side = 1, at = xvec, labels = xlab)
    errbar(x = xvec, y = power.rnd.mean, yplus = power.rnd.plus, yminus = power.rnd.minus
           , xaxt = "n", xlab = "caliper size", ylab = "power"
           , cex = cex.mult, cex.axis = cex.mult, cex.lab = cex.mult, cex.main = cex.mult
           , lty = 2, type = "b"
           , add = TRUE
           )
    legend("bottomright", legend = c("matched subsample", "random subsample"), lty = 1:2)
  }
  
}

mlr <- function(tr, Z.i = NULL, Z.o = mlr.generate.Z.o(Z.i), psm = TRUE
  , caliper.vec = c(0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 5.0, Inf)
  , ...) {
  
  Z.match <- cbind(Z.i, Z.o) # perhaps we can apply svd to Z.match before feeding it to logistic reg
  idx.list <- lapply(caliper.vec, function(caliper) {
    mlr.match(tr = tr, X = Z.match, psm = psm, verbose = FALSE, caliper = caliper, ...)
  })

  delta.x.1 <- 3.0
  delta.x.2 <- 1 * delta.x.1
  xvec.orig <- caliper.vec
  xvec <- xvec.orig
  xvec[is.infinite(xvec)] <- max(xvec[!is.infinite(xvec)]) * delta.x.1
  xvec <- c(xvec, max(xvec) * delta.x.2) # for no matching
  xlab <- format(xvec.orig, digits = 2); xlab[is.infinite(xvec.orig)] <- "no caliper"; xlab <- c(xlab, "no match")
  xtitle <- "caliper size" # this should also come from mlr
  
  ret <- list(tr = tr, Z.i = Z.i, Z.o = Z.o, idx.list = idx.list
              #, caliper.vec = caliper.vec
              , xvec = xvec, xlab = xlab, xtitle = xtitle
              )
  class(ret) <- "mlr"
  
  return (ret)
}

mlr.combine.bias.variance <- function(tr, bvmat, orsq.min = 1e-3, orsq.max = 1e+0, n.orsq = 100) {
  if (length(bvmat) == NROW(bvmat)) bvmat <- t(as.matrix(bvmat))
  else bvmat <- as.matrix(bvmat)
  n.method <- nrow(bvmat)
  
  nt <- length(which(tr == 1))
  nc <- length(which(tr == 0))
  var.min <- 1/nc + 1/nt
  
  orsq.seq <- exp(seq(from = log(orsq.min), to = log(orsq.max), length.out = n.orsq))
  errmat <- array(NA, dim = c(n.orsq, n.method))
  biassq.mat <- array(NA, dim = c(n.orsq, n.method))
  which.min.vec <- rep(NA, n.method)
  for (i in 1:n.orsq) {
    biassq.mat[i, ] <- orsq.seq[i] * bvmat[, 1]^2
    errmat[i, ] <- c(orsq.seq[i] * bvmat[, 1]^2 + bvmat[, 2])
    which.min.vec[i] <- which.min(errmat[i, ])
  }
  
  ret <- list(orsq.vec = orsq.seq, errmat = errmat, biassq.mat = biassq.mat, which.min.vec = which.min.vec)
  
  return (ret)
}

