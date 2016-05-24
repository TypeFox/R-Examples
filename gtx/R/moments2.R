is.moments2 <- function(object) {
  if (class(object) != "moments2") return (FALSE)
  if (!is.matrix(object)) return (FALSE)
  if (nrow(object) != ncol(object)) return (FALSE)
  if (any(rownames(object) != colnames(object))) return (FALSE)
  return (TRUE)
}

make.moments2 <- function(params, phenolist, snpdata, weightvar = NULL) {
  stopifnot(is.data.frame(params))
  stopifnot(all(c("snp", "coded.allele") %in% names(params)))
  stopifnot(is.snpdata(snpdata))
  if (! "ONE" %in% names(snpdata$data)) snpdata$data$ONE <- 1 # ensures no column called "ONE" would be overwritten
  stopifnot(all(snpdata$data$ONE == 1))
  stopifnot(all(phenolist %in% names(snpdata$data)))
  if (!is.null(weightvar)) stopifnot(length(weightvar) == 1)
  if (!is.null(weightvar)) stopifnot(weightvar %in% names(snpdata$data))	
  snpdata <- align.snpdata.coding(params, snpdata, missing.snp = "okay")$snpdata
  tmp <- subset(snpdata$data, select = c("ONE", paste(params$snp, params$coded.allele, sep = "_"), phenolist, weightvar))
  if (!all(sapply(names(tmp), function(nn) is.numeric(tmp[[nn]])))) {
    stop("non-numeric data encountered")
  }
  if (!is.null(weightvar)) widx <- match(weightvar, names(tmp))
  x <- as.matrix(subset(tmp, !apply(is.na(tmp), 1, any)))
  ## note we extract columns of interest first, THEN hard drops any rows with missing data
  if (is.null(weightvar)) {
    xtx <- t(x) %*% x
    class(xtx) <- "moments2"
    attr(xtx, "vscale") <- NULL
    return(xtx)
  } else {
    xtwx <- t(apply(x, 2, function(xcol) xcol*x[ , widx])) %*% x
    class(xtwx) <- "moments2"
    attr(xtwx, "vscale") <- 1
    return(xtwx)
  }
}

lm.moments2 <- function(xtx, leftvar, rightvars, n = NULL) {
  return(est.moments2(xtx, leftvar, rightvars, n, vscale = NULL))
}

est.moments2 <- function(xtwx, leftvar, rightvars, n = NULL, vscale = NULL) {
  ## estimating vscale from data assumes normal/identity model
  ## if using vscale!=1 check definition vis-a-vis standard GLM notation, maybe inverse?

  if(is.moments2(xtwx)) {
    vscale.arg <- vscale
    vscale <- attr(xtwx, "vscale")
    if ((is.null(vscale.arg) && !is.null(vscale)) ||
        (!is.null(vscale.arg) && is.null(vscale)) ||
        (!is.null(vscale.arg) && !is.null(vscale) && vscale.arg != vscale)) {
      warning("using vscale=", if (is.null(vscale)) "NULL" else vscale, " from xt(w)x attribute", 
              "instead of vscale=" , if (is.null(vscale.arg)) "NULL" else vscale.arg, " from function argument")
    }
  } else {
    ## these checks only needed if !is.moments2(xtwx)
    stopifnot(is.matrix(xtwx))
    stopifnot(nrow(xtwx) == ncol(xtwx))
    stopifnot(all(rownames(xtwx) == colnames(xtwx)))
  }
  if (is.null(n)) stopifnot("ONE" %in% rownames(xtwx))
  stopifnot(length(leftvar) == 1)
  stopifnot(leftvar %in% rownames(xtwx))
  stopifnot(length(rightvars) >= 1)
  stopifnot(all(rightvars %in% rownames(xtwx)))

  ## n inferred this way is only right for NLM
  if (is.null(vscale) && is.null(n)) n <- xtwx["ONE", "ONE"]
  p <- length(rightvars)
  stopifnot(n - p >= 1) # residual degrees of freedom; nothing happens if n==NULL

  lidx <- match(leftvar, rownames(xtwx))
  
  ridx <- match(rightvars, rownames(xtwx))
  myxtwxi <- try(solve(xtwx[ridx, ridx, drop = FALSE]), silent = TRUE)
  if (class(myxtwxi) == "matrix") { # simple case X'X is nonsingular
    myxtwy <- xtwx[ridx, lidx, drop = FALSE]
    betahat <- drop(myxtwxi %*% myxtwy)
    if (is.null(vscale)) {
      ssr <- as.double(xtwx[lidx, lidx, drop = TRUE] - t(myxtwy) %*% myxtwxi %*% myxtwy)
      vscale <- ssr/(n - p)
    } else {
      ssr <- NA
    }
    se <- sqrt(diag(myxtwxi)*vscale)
    prec <- xtwx[ridx, ridx, drop = FALSE]/vscale
    return(list(betahat = betahat, se = se, prec = prec, ssr = ssr))
  } else if (class(myxtwxi) == "try-error") {
    stopifnot(length(ridx) > 1)     #### if singular with one variable then only coefficient is nonidentifiable
    uridx <- ridx[1] # useable right indices
    for (widx in ridx[-1]) {
      if (class(try(solve(xtwx[c(uridx, widx), c(uridx, widx), drop = FALSE]), silent = TRUE)) == "matrix"){
        uridx <- c(uridx, widx)
##      } else {
##        warning("dropping ", rownames(xtwx)[widx], " due to singularity")
      }
    }
    ## warn if setdiff(ridx, uridx) is non-empty
    ndrop <- length(setdiff(ridx, uridx))
    if (ndrop > 0) warning (ndrop, " variables with non-identifiable coefficients dropped")
    ## something odd happened if ndrop==0!
    ## NA fill vectors & matrices of right size
    betahat <- rep(NA, length(rightvars)) ## or zero if we take limit with prior
    se <- rep(NA, length(rightvars)) ## or Inf if we take limit with prior
    #prec <- matrix(rep(NA, length(rightvars)^2), nrow = length(rightvars))

    myxtwxi <- solve(xtwx[uridx, uridx, drop = FALSE])
    myxtwy <- xtwx[uridx, lidx, drop = FALSE]
    betahat[match(uridx, ridx)] <- myxtwxi %*% myxtwy
    if (is.null(vscale)) {
      ssr <- as.double(xtwx[lidx, lidx, drop = TRUE] - t(myxtwy) %*% myxtwxi %*% myxtwy)
      vscale <- ssr/(n - p)
    } else {
      ssr <- NA
    }
    se[match(uridx, ridx)] <- sqrt(diag(myxtwxi)*vscale)
    ##prec[match(uridx, ridx), match(uridx, ridx)] <- xtwx[uridx, uridx, drop = FALSE]/vscale
    prec <- xtwx[ridx, ridx, drop = FALSE]/vscale # IS THIS RIGHT?
    return(list(betahat = betahat, se = se, prec = prec, ssr = ssr))
  } else {
    stop("internal error inverting X'X")
  }
}
       
combine.moments2 <- function(xtxlist, fixed) {
  m <- length(xtxlist)
  stopifnot(!("ONE" %in% fixed))
  stopifnot(length(intersect(names(xtxlist), fixed)) == 0)
  stopifnot(all(sapply(names(xtxlist), function(ii) is.matrix(xtxlist[[ii]]))))
  stopifnot(all(sapply(names(xtxlist), function(ii) all(rownames(xtxlist[[ii]]) == colnames(xtxlist[[ii]])))))
  stopifnot(all(sapply(names(xtxlist), function(ii) all(c("ONE", fixed) %in% rownames(xtxlist[[ii]])))))
  
  free <- list()
  for (ii in names(xtxlist)) free[[ii]] <- setdiff(rownames(xtxlist[[ii]]), c("ONE", fixed))

  si <- sapply(1:m, function(ii) length(free[[ii]]))
  stopifnot(length(unlist(free)) == sum(si))
  osi <- m + length(fixed) + c(0, cumsum(si)[-m]) # offset-1 in xtx where covariates for i-th study start

  xtx <- matrix(rep(0, (m + length(fixed) + sum(si))^2), ncol = m + length(fixed) + sum(si))
  rownames(xtx) <- 1:(m + length(fixed) + sum(si))
  rownames(xtx)[1:m] <- paste(names(xtxlist), "ONE", sep = "_")
  rownames(xtx)[m + 1:length(fixed)] <- fixed
  for (ii in which(si > 0)) rownames(xtx)[osi[ii] + 1:si[ii]] <- paste(names(xtxlist)[ii], free[[ii]], sep = "_")
  colnames(xtx) <- rownames(xtx)
    
  for (ii in 1:m) xtx[ii, ii] <- xtxlist[[ii]]["ONE", "ONE", drop = FALSE]
  for (ii in 1:m) xtx[ii, m + 1:length(fixed)] <- xtxlist[[ii]]["ONE", fixed, drop = FALSE]
  for (ii in 1:m) xtx[m + 1:length(fixed), ii] <- xtxlist[[ii]][fixed, "ONE", drop = FALSE]
  for (ii in 1:m) xtx[m + 1:length(fixed), m + 1:length(fixed)] <- xtx[m + 1:length(fixed), m + 1:length(fixed), drop = FALSE] + xtxlist[[ii]][fixed, fixed, drop = FALSE]
  
  for (ii in which(si > 0)) {
    xtx[ii, osi[ii] + 1:si[ii]] <- xtxlist[[ii]]["ONE", free[[ii]], drop = FALSE]
    xtx[osi[ii] + 1:si[ii], ii] <- xtxlist[[ii]][free[[ii]], "ONE", drop = FALSE]
    xtx[m + 1:length(fixed), osi[ii] + 1:si[ii]] <- xtxlist[[ii]][fixed, free[[ii]], drop = FALSE]
    xtx[osi[ii] + 1:si[ii], m + 1:length(fixed)] <- xtxlist[[ii]][free[[ii]], fixed, drop = FALSE]
    xtx[osi[ii] + 1:si[ii], osi[ii] + 1:si[ii]] <- xtxlist[[ii]][free[[ii]], free[[ii]], drop = FALSE]
  }
  return(xtx)
}
