fmri.smooth <- function (spm, hmax = 4, adaptation = "aws",
    lkern = "Gaussian", skern = "Plateau", weighted = TRUE, ...)
{
    cat("fmri.smooth: entering function\n")
    args <- sys.call()
    args <- c(spm$call,args)
    if (!(tolower(adaptation) %in% c("none", "aws", "segment"))) {
        adaptation <- "aws"
    }
    ladjust <- if ("ladjust" %in% names(list(...)))
        list(...)[["ladjust"]]
    else 1
    delta <- if ("delta" %in% names(list(...)))
        list(...)[["delta"]]
    else 0
    alpha <- if ("alpha" %in% names(list(...)))
        list(...)[["alpha"]]
    else 0.05
    propagation <- if ("propagation" %in% names(list(...)))
        list(...)[["propagation"]]
    else FALSE
    restricted <- if ("restricted" %in% names(list(...)))
        as.logical(list(...)[["restricted"]])
    else FALSE
    if (!("fmrispm" %in% class(spm))) {
        warning("fmri.smooth: data not of class <fmrispm>. Try to proceed but strange things may happen")
    }
    if (!is.null(attr(spm, "smooth"))) {
        warning("fmri.smooth: Parametric Map seems to be smoothed already!")
    }
    variance <- spm$var
    variance[variance == 0] <- 1e+20
    if (is.null(spm$weights)) {
        weights <- c(1, 1, 1)
    }
    else {
        weights <- spm$weights
    }
    if (is.null(spm$bw)) {
        bw <- rep(0, 3)
    }
    else {
        bw <- spm$bw
    }
    cat("fmri.smooth: smoothing the Statistical Parametric Map\n")
    ttthat <- switch(tolower(adaptation), 
                     aws = vaws3D(y = spm$cbeta,
                                  sigma2 = variance, 
                                  hmax = hmax, 
                                  mask = spm$mask, 
                                  wghts = weights,
                                  h0 = bw, 
                                  vwghts = spm$vwghts, 
                                  lkern = lkern, 
                                  skern = skern,
                                  weighted = weighted, 
                                  res = spm$res, 
                                  resscale = spm$resscale,
                                  ddim = spm$dim, 
#                                  df = spm$df,
                                  ladjust = ladjust, 
                                  testprop = propagation),
                     fullaws = vaws3Dfull(y = spm$cbeta, 
                                          sigma2 = variance,
                                          hmax = hmax, 
                                          mask = spm$mask, 
                                          wghts = weights, 
                                          vwghts = spm$vwghts,
                                          lkern = lkern, 
                                          skern = skern, 
                                          weighted = weighted,
                                          res = spm$res, 
                                          resscale = spm$resscale, 
                                          ddim = spm$dim,
#                                          df = spm$df,
                                          ladjust = ladjust, 
                                          testprop = propagation),
                     none = vaws3D(y = spm$cbeta,
                                   sigma2 = variance, 
                                   hmax = hmax, 
                                   mask = spm$mask,
                                   qlambda = 1, 
                                   wghts = weights, 
                                   h0 = bw, 
                                   vwghts = spm$vwghts,
                                   lkern = lkern, 
                                   skern = skern, 
                                   weighted = weighted,
                                   res = spm$res, 
                                   resscale = spm$resscale, 
#                                   df = spm$df,
                                   ddim = spm$dim,
                                   ladjust = ladjust),
                     segment = segm3D(y = spm$cbeta,
                                      sigma2 = variance, 
                                      hmax = hmax, 
                                      mask = spm$mask,
                                      wghts = weights, 
                                      df = spm$df, 
                                      h0 = bw, 
                                      weighted = weighted,
                                      res = spm$res, 
                                      resscale = spm$resscale, 
                                      ddim = spm$dim,
                                      ladjust = ladjust, 
                                      delta = delta, 
                                      alpha = alpha,
                                      restricted = restricted))
    cat("\n")
    cat("fmri.smooth: determine local smoothness\n")
    if (is.null(ttthat$scorr)) {
        bw <- get3Dh.gauss(ttthat$vred, weights)
    }
    else {
        bw <- optim(c(2, 2, 2), corrrisk, method = "L-BFGS-B",
            lower = c(0.25, 0.25, 0.25), upper = c(6, 6, 6),
            lag = c(5, 5, 3), data = ttthat$scorr)$par
        bw[bw <= 0.25] <- 0
        dim(bw) <- c(1, 3)
    }
    rxyz <- c(resel(1, bw[, 1]), resel(1, bw[, 2]), resel(1,
        bw[, 3]))
    dim(rxyz) <- c(dim(bw)[1], 3)
    bw0 <- get3Dh.gauss(ttthat$vred0, weights)
    rxyz0 <- c(resel(1, bw0[, 1]), resel(1, bw0[, 2]), resel(1,
        bw0[, 3]))
    dim(rxyz0) <- c(dim(bw0)[1], 3)
    cat("fmri.smooth: exiting function\n")
    if (length(dim(ttthat$theta)) == 3)
        dim(ttthat$theta) <- c(dim(ttthat$theta), 1)
    if (dim(ttthat$theta)[4] == 1) {
        z <- list(cbeta = ttthat$theta[, , , 1], var = ttthat$var,
            rxyz = rxyz, rxyz0 = rxyz0, scorr = spm$scorr, weights = spm$weights,
            vwghts = spm$vwghts, bw = bw, hmax = ttthat$hmax,
            dim = spm$dim, hrf = spm$hrf, segm = ttthat$segm,
            mask = ttthat$mask, call = args)
    }
    else {
        z <- list(cbeta = ttthat$theta, var = ttthat$var, rxyz = rxyz,
            rxyz0 = rxyz0, scorr = spm$scorr, weights = spm$weights,
            vwghts = spm$vwghts, bw = bw, hmax = ttthat$hmax,
            dim = spm$dim, hrf = spm$hrf, segm = ttthat$segm,
            mask = ttthat$mask, call = args)
    }
    if (adaptation == "segment") {
      class(z) <- c( "fmrisegment", "fmridata")
      z$alpha <- alpha
      z$delta <- delta
    } else {
      class(z) <- c("fmridata", "fmrispm")
    }
    z$roixa <- spm$roixa
    z$roixe <- spm$roixe
    z$roiya <- spm$roiya
    z$roiye <- spm$roiye
    z$roiza <- spm$roiza
    z$roize <- spm$roize
    z$roit <- spm$roit
    z$header <- spm$header
    z$format <- spm$format
    z$dim0 <- spm$dim0
    z$scorr <- ttthat$scorr
    z$call <- args
    attr(z, "file") <- attr(spm, "file")
    attr(z, "white") <- attr(spm, "white")
    attr(z, "design") <- attr(spm, "design")
    attr(z, "residuals") <- attr(spm, "residuals")
    if (!is.null(attr(spm, "smooth"))) {
        attr(z, "smooth") <- paste("Already smoothed before:\n",
            attr(spm, "smooth"), "\nnow with:\n  adaptation  :",
            adaptation, "\n  bandwidth :", signif(hmax,
                3), "\n  lkern     :", lkern, "\n  skern     :",
            skern, "\n")
    }
    else {
        attr(z, "smooth") <- paste("Smoothed with:\n  adaptation  :",
            adaptation, "\n  bandwidth :", signif(hmax,
                3), "\n  lkern     :", lkern, "\n  skern     :",
            skern, "\n")
    }
    invisible(z)
}

fmri.pvalue <- function(spm, mode="basic", delta=NULL, na.rm=FALSE, minimum.signal=0, alpha=0.05 ) {
    args <- sys.call()
    args <- c(spm$call,args)
  cat("fmri.pvalue: entering function\n")

  if (!("fmrispm" %in% class(spm)) ) {
    warning("fmri.pvalue: data not of class <fmrispm>. Try to proceed but strange things may happen")
  }

  if (!is.null(attr(spm, "smooth"))) {
    if (!is.null(attr(spm, "residuals"))) {
      type <- "t"
      df <- spm$df 
      if(is.null(df)) df <- abs(diff(dim(attr(spm, "design"))))
    } else {
      type <- "norm"
      df <- 1000 # this is actually not needed, placeholder
    }
  } else {
    type <- "t"
    df <- spm$df
  }

  if (length(dim(spm$cbeta)) < 4) {

    stat <- (spm$cbeta-minimum.signal)/sqrt(spm$var)
    dim(stat) <- prod(spm$dim[1:3])
    cat("fmri.pvalue: calculate treshold and p-value method:",mode,"\n")
    if (mode == "local") {
      thresh <- threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],spm$rxyz[,1],spm$rxyz[,2],spm$rxyz[,3],type=type,df=df)
      pv <- pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],spm$rxyz[,1],spm$rxyz[,2],spm$rxyz[,3],type=type,df=df)
    } else if (mode == "global") {
      rxyz <- c(median(spm$rxyz[,1]),median(spm$rxyz[,2]),median(spm$rxyz[,3]))
      thresh <- threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],rxyz[1],rxyz[2],rxyz[3],type=type,df=df)
      pv <- pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],rxyz[1],rxyz[2],rxyz[3],type=type,df=df)
     } else if (mode == "FDR") {
      pv <- 1-switch(type,"norm"=pnorm(stat),"t"=pt(stat,df))
      dim(pv) <- spm$dim[1:3]
      ind <- fdr(pv,alpha)
      thresh <- min(stat[ind])
   } else {
      if ("rxyz0" %in% names(spm)) {
        rxyz0 <- c(median(spm$rxyz0[,1]),median(spm$rxyz0[,2]),median(spm$rxyz0[,3]))
      } else {
        rxyz0 <- c(median(spm$rxyz[,1]),median(spm$rxyz[,2]),median(spm$rxyz[,3]))
      }        
      thresh <- threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],rxyz0[1],rxyz0[2],rxyz0[3],type=type,df=df)
      pv <- pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],rxyz0[1],rxyz0[2],rxyz0[3],type=type,df=df)
    }

  } else if (!is.null(delta)) {

    l1 <- sqrt(spm$vwghts[2]/spm$vwghts[1]) * delta[1]
    l2 <- sqrt(spm$vwghts[2]/spm$vwghts[1]) * delta[2]
    theta1 <- atan(l1)
    theta2 <- atan(l2)
    t1 <- spm$cbeta[,,,1]/sqrt(spm$var * spm$vwghts[1])
    t2 <- spm$cbeta[,,,2]/sqrt(spm$var * spm$vwghts[2])
    ratio <- t2/t1
    ratio[t1==0] <- l2 + 1
    w1 <- (t1 + t2 * l1) / sqrt(1+l1^2)
    w2 <- (t1 + t2 * l2) / sqrt(1+l2^2)
    w3 <- (t1 > 0) * (l1 <= ratio) * (ratio <= l2) * sqrt(t1^2 + t2^2)
    stat <- pmax(w1,w2,w3)
    dim(stat) <- prod(spm$dim[1:3])
    cat("fmri.pvalue: calculate p-value method:",mode,"\n")
    if (mode == "local") {
      thresh <-
        threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],spm$rxyz[,1],spm$rxyz[,2],spm$rxyz[,3],type="norm",cone=theta2-theta1)
      pv <-
        pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],spm$rxyz[,1],spm$rxyz[,2],spm$rxyz[,3],type="norm",cone=theta2-theta1)
    } else if (mode == "global") {
      rxyz <- c(median(spm$rxyz[,1]),median(spm$rxyz[,2]),median(spm$rxyz[,3]))
      thresh <-
        threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],rxyz[1],rxyz[2],rxyz[3],type="norm",cone=theta2-theta1)
      pv <-
        pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],rxyz[1],rxyz[2],rxyz[3],type="norm",cone=theta2-theta1)
    } else {
      if ("rxyz0" %in% names(spm)) {
        rxyz0 <- c(median(spm$rxyz0[,1]),median(spm$rxyz0[,2]),median(spm$rxyz0[,3]))
      } else {
        rxyz0 <- c(median(spm$rxyz[,1]),median(spm$rxyz[,2]),median(spm$rxyz[,3]))
      }
      thresh <-
        threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],rxyz0[1],rxyz0[2],rxyz0[3],type="norm",cone=theta2-theta1)
      pv <-
        pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],rxyz0[1],rxyz0[2],rxyz0[3],type="norm",cone=theta2-theta1)
    }

  } else {

    stat <- spm$cbeta[,,,1]^2/spm$var + spm$cbeta[,,,2]^2/spm$var/spm$vwghts[2]  # Wert der Statistik
    dim(stat) <- prod(spm$dim[1:3])
    cat("fmri.pvalue: calculate treshold and p-value method:",mode,"\n")
    if (mode == "local") {
      thresh <-
        threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],spm$rxyz[,1],spm$rxyz[,2],spm$rxyz[,3],type="chisq",df=2)
      pv <-
        pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],spm$rxyz[,1],spm$rxyz[,2],spm$rxyz[,3],type="chisq",df=2)
    } else if (mode == "global") {
      rxyz <- c(median(spm$rxyz[,1]),median(spm$rxyz[,2]),median(spm$rxyz[,3]))
      thresh <-
        threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],rxyz[1],rxyz[2],rxyz[3],type="chisq",df=2)
      pv <-
        pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],rxyz[1],rxyz[2],rxyz[3],type="chisq",df=2)
    } else {
      if ("rxyz0" %in% names(spm)) {
        rxyz0 <- c(median(spm$rxyz0[,1]),median(spm$rxyz0[,2]),median(spm$rxyz0[,3]))
      } else {
        rxyz0 <- c(median(spm$rxyz[,1]),median(spm$rxyz[,2]),median(spm$rxyz[,3]))
      }
      thresh <-
        threshold(0.2,spm$dim[1],spm$dim[2],spm$dim[3],rxyz0[1],rxyz0[2],rxyz0[3],type="chisq",df=2)
      pv <-
        pvalue(stat,spm$dim[1],spm$dim[2],spm$dim[3],rxyz0[1],rxyz0[2],rxyz0[3],type="chisq",df=2)
    }

  }
  cat("fmri.pvalue: thresholding\n")
  mask <- rep(TRUE,length=prod(spm$dim[1:3]))
  mask[stat < thresh] <- FALSE
  mask <- mask&spm$mask
  pv[!mask] <- 1
  dim(pv) <- spm$dim[1:3]

  if (na.rm) {
    pv[spm$var > 9e19] <- 1
  }
  
  cat("fmri.pvalue: exiting function\n")

  z <- list(pvalue = pv, weights = spm$weights, dim = spm$dim, hrf = spm$hrf)
  
  class(z) <- c("fmridata","fmripvalue")

  z$roixa <- spm$roixa
  z$roixe <- spm$roixe
  z$roiya <- spm$roiya
  z$roiye <- spm$roiye
  z$roiza <- spm$roiza
  z$roize <- spm$roize
  z$roit <- spm$roit
  z$header <- spm$header
  z$format <- spm$format
  z$dim0 <- spm$dim0
  z$args <- z$args

  attr(z, "file") <- attr(spm, "file")
  attr(z, "white") <- attr(spm, "white")
  attr(z, "design") <- attr(spm, "design")
  if (is.null(attr(spm, "smooth"))) {
    attr(z, "smooth") <- "Not smoothed"
  } else {
    attr(z, "smooth") <- attr(spm, "smooth")
  }
  attr(z, "mode") <- paste("Threshold mode:",mode,"\n")
   
  z
}

