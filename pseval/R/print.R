
#' Plot summary statistics for a psdesign object
#'
#' Plot the vaccine efficacy or another contrast of risk versus S.1 for an
#' estimated psdesign object
#'
#' @param x A \link{psdesign} object that contains a risk model, integration
#'   model, and valid estimates
#' @param t For time to event outcomes, a fixed time \code{t} may be provided to
#'   compute the cumulative distribution function. If not, the restricted mean
#'   survival time is used. Omit for binary outcomes.
#' @param contrast Name of contrast function to plot. \code{"VE"} for vaccine
#'   efficacy = 1 - risk_1(s)/risk_0(s), \code{"RR"} for relative risk =
#'   risk_1(s)/risk_0(s), \code{"logRR"} for log of the relative risk,
#'   \code{"risk"} for the risk in each treatment arm, and \code{"RD"} for the
#'   risk difference = risk_1(s) - risk_0(s). You can also pass a custom
#'   function directly as long as it takes 2 vectors as input (risk0 and risk1)
#'   and returns 1 vector of the same length.
#' @param sig.level Significance level used for confidence bands on the VE
#'   curve. This is only used if bootstrapped estimates are available.
#' @param CI.type Character string, "pointwise" for pointwise confidence
#'   intervals, and "band" for simultaneous confidence band.
#' @param n.samps Number of samples to use over the range of S.1 for plotting
#'   the curve
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param col Vector of integers specifying colors for each curve.
#' @param lty Vector of integers specifying linetypes for each curve.
#' @param lwd Vector of numeric values for line widths.
#' @param ... Other arguments passed to \link{plot}
#'
#' @export

plot.psdesign <- function(x, t, contrast = "VE", sig.level = .05, CI.type = "band", n.samps = 500, xlab = "S.1", ylab = contrast, col = 1, lty = 1, lwd = 1, ...){


  if(contrast == "risk"){
    contrast0 <- "VE"
  } else{
    contrast0 <- contrast
  }
  VE.me <- calc_risk(x, contrast0, t, sig.level = sig.level, CI.type = CI.type, n.samps = n.samps)
  n.curve.base <- ifelse(contrast == "risk", 2, 1)
  ncurve <- ifelse("Y.boot.se" %in% colnames(VE.me), n.curve.base * 3, n.curve.base)

  ## some logic taken from plot.survfit


    if (length(lty)==1 && is.numeric(lty))
      lty <- rep(c(lty, lty+1, lty+1), 2)
    else if (length(lty) < ncurve)
      lty <- rep(rep(lty, each=3), length.out=(ncurve*3))
    else lty <- rep(lty, length.out= ncurve*3)

    if (length(col) <= ncurve) col <- rep(rep(col, each=3), length.out=3*ncurve)
    else col <- rep(col, length.out=3*ncurve)

    if (length(lwd) <= ncurve) lwd <- rep(rep(lwd, each=3), length.out=3*ncurve)
    else lwd <- rep(lwd, length.out=3*ncurve)


  mainme <- switch(contrast,
                   risk = parse(text = "list(R1, R0)"))


  if(is.factor(VE.me[, 1])){
    envir <- unique(VE.me)
  } else envir <- VE.me

  if(contrast == "risk"){

    rlist <- eval(mainme, envir = envir)
    plot(rlist[[1]] ~ envir[, 1], type = 'l', col = col[1], lty = lty[1], lwd = lwd[1],
         ylab = ylab, xlab = xlab, ...)
    if(is.factor(VE.me[, 1])){

      segments(rep.int(as.integer(envir[, 1]), 2) - .4, rep.int(rlist[[2]], 2),
               x1 = rep.int(as.integer(envir[, 1]), 2) + .4, col = col[4], lty = lty[4], lwd = lwd[4])

    } else {

      lines(rlist[[2]] ~ envir[, 1], type = 'l', col = col[4], lty = lty[4], lwd = lwd[4])

    }

  } else {

    plot(envir$Y ~ envir$S.1, col = col[1], lty = lty[1],
         lwd = lwd[1], type = 'l', ylab = ylab, xlab = xlab, ...)

  }


  if("bootstraps" %in% names(x)){

    lnme <- switch(contrast,
                   risk = parse(text = paste("list(", paste(sapply(c("R1.lower.", "R0.lower."),
                                                                   function(x) grep(x, colnames(VE.me), fixed = TRUE, value = TRUE)), collapse = ", "), ")")))

    unme <- switch(contrast,
                   risk = parse(text = paste("list(", paste(sapply(c("R1.upper.", "R0.upper."),
                                                                   function(x) grep(x, colnames(VE.me), fixed = TRUE, value = TRUE)), collapse = ", "), ")")))

    if(contrast == "risk"){

      rlistl <- eval(lnme, envir = envir)
      rlistu <- eval(unme, envir = envir)

      if(is.factor(VE.me[, 1])){

        segments(rep.int(as.integer(envir[, 1]), 2) - .4, rep.int(rlistu[[1]], 2),
                 x1 = rep.int(as.integer(envir[, 1]), 2) + .4, col = col[2], lty = lty[2], lwd = lwd[2], ...)
        segments(rep.int(as.integer(envir[, 1]), 2) - .4, rep.int(rlistu[[2]], 2),
                 x1 = rep.int(as.integer(envir[, 1]), 2) + .4, col = col[5], lty = lty[5], lwd = lwd[5], ...)
        segments(rep.int(as.integer(envir[, 1]), 2) - .4, rep.int(rlistl[[1]], 2),
                 x1 = rep.int(as.integer(envir[, 1]), 2) + .4, col = col[3], lty = lty[3], lwd = lwd[3], ...)
        segments(rep.int(as.integer(envir[, 1]), 2) - .4, rep.int(rlistl[[2]], 2),
                 x1 = rep.int(as.integer(envir[, 1]), 2) + .4, col = col[6], lty = lty[6], lwd = lwd[6], ...)

      } else {

        lines(rlistu[[1]] ~ envir[, 1], type = 'l', col = col[2], lty = lty[2], lwd = lwd[2])
        lines(rlistu[[2]] ~ envir[, 1], type = 'l', col = col[5], lty = lty[5], lwd = lwd[5])
        lines(rlistl[[1]] ~ envir[, 1], type = 'l', col = col[3], lty = lty[3], lwd = lwd[3])
        lines(rlistl[[2]] ~ envir[, 1], type = 'l', col = col[6], lty = lty[6], lwd = lwd[6])

      }

    } else {

      if(is.factor(VE.me[, 1])){

        segments(rep.int(as.integer(envir[, 1]), 2) - .4, rep.int(envir[, 6], 2),
                 x1 = rep.int(as.integer(envir[, 1]), 2) + .4, col = col[2], lty = lty[2], lwd = lwd[2], ...)
        segments(rep.int(as.integer(envir[, 1]), 2) - .4, rep.int(envir[, 7], 2),
                 x1 = rep.int(as.integer(envir[, 1]), 2) + .4, col = col[3], lty = lty[3], lwd = lwd[3], ...)

      } else {

        lines(envir[, 6] ~ envir[, 1], col = col[2], lty = lty[2], lwd = lwd[2], type = 'l')
        lines(envir[, 7] ~ envir[, 1], col = col[3], lty = lty[3], lwd = lwd[3], type = 'l')
      }
    }


    }


}


#' Concisely print information about a psdesign object
#'
#' @param x An object of class \link{psdesign}
#' @param digits Number of significant digits to display
#' @param sig.level Significance level to use for computing bootstrapped confidence intervals
#' @param ... Currently unused
#' @export
#'
print.psdesign <- function(x, digits = 3, sig.level = .05, ...){

  objs <- names(x)
  pout <- NULL

  cat("Augmented data frame: ", nrow(x$augdata), " obs. by ", ncol(x$augdata), " variables. \n")
  print(head(x$augdata), digits = digits)

  cat("\nEmpirical VE: ", round(empirical_VE(x), digits), "\n")

  cat("\nMapped variables: \n")
  temp <- lapply(names(x$mapping), function(ln){

    lnm <- x$mapping[[ln]]

    if(length(lnm) == 1){
      cat("\t", ln, " -> ", lnm, "\n")
    } else if(lnm[1] == "Surv") {

      cat("\t", ln, " -> ", paste0(lnm[1], "(", lnm[2], ", ", lnm[3], ")"), "\n")
    } else cat("\t", ln, " -> ", paste0(lnm, collapse = ", "), "\n")

  })

  cat("\nIntegration models: \n")
  if(!"integration.models" %in% objs) {

    cat("\tNone present, see ?add_integration for information on integration models.\n")

  } else {

    for(j in names(x$integration.models)){
      cat("\t integration model for ", j, ":\n")
      tyj <- x$integration.models[[j]]$model
      cat("\t\t", paste0("integrate_", tyj$model, "("))
      cat(paste(sapply(names(tyj$args[-1]), function(nj) paste0(nj, " = ", tyj$args[-1][nj])), collapse = ", "), ")\n")
    }

  }

  cat("\nRisk models: \n")
  if(!"risk.model" %in% objs) {
    cat("\tNone present, see ?add_riskmodel for information on risk models.\n")
  } else {

    tyj <- x$risk.model
    cat("\t", paste0("risk_", tyj$model, "("))
    cat(paste(sapply(names(tyj$args[-1]), function(nj) paste0(nj, " = ", tyj$args[-1][nj])), collapse = ", "), ")\n\n")

  }

  if(!"estimates" %in% objs){

    cat("\tNo estimates present, see ?ps_estimate.\n")
  } else {
    cat("Estimated parameters:\n")
    print(x$estimates$par, digits = digits)
    cat("\tConvergence: ", x$estimates$convergence == 0, "\n\n")

  }

  if(!"bootstraps" %in% objs){

    cat("\tNo bootstraps present, see ?ps_bootstrap.\n")
  } else {

    cat("Bootstrap replicates:\n")
    sbs <- summarize_bs(x$bootstraps, x$estimates$par, sig.level = sig.level, CI.type = "pointwise")
    print(sbs$table, digits = digits)
    cat("\n\t Out of", sbs$conv[1], "bootstraps, ", sbs$conv[2], "converged (", round(100 * sbs$conv[2]/sbs$conv[1], 1), "%)\n")

    ## WEM test
    wem <- wem_test(x)
    cat("\n\t Test for wide effect modification on", wem$df, ifelse(wem$df > 1, "degrees", "degree"),
        "of freedom. 2-sided p value ",
        ifelse(wem$p.value < .0001, "< .0001", paste0("= ", round(wem$p.value, 4))), "\n")


    pout$wem.test <- wem
    pout$boot.table <- sbs

  }

  invisible(pout)


}

#' Summary method for psdesign objects
#'
#' @return Invisibly returns the printed table, along with the three estimates of vaccine efficacy. The empirical VE is 1 minus the relative risk comparing the treatment arm to the control arm. The risk is estimated as the proportion in the binary outcome case, or with the Kaplan-Meier estimate at the restricted mean survival in the time-to-event case. The marginal VE estimate is the VE estimate under the specified parametric risk model, ignoring the effect of S.1. The model based average VE is the VE estimate from the specified risk model, averaged over the distribution of S.1. The point of displaying these three is to assess the validity of the parametric model, and to assess the validity of the model estimation. Wild differences among these estimates may indicate problems with the model or convergence.
#'
#' @param object An object of class \link{psdesign}
#' @param digits Number of significant digits to display
#' @param sig.level Significance level to use for computing bootstrapped confidence intervals
#' @param ... Currently unused
#'
#' @export
#'
summary.psdesign <- function(object, digits = 3, sig.level = .05, ...){

  pout <- print(object, digits = digits, sig.level = sig.level)

  ## compute marginal model and summarize VE

  if("risk.model" %in% names(object)){
    pdat <- object$risk.model$args
    pdat$model <- Y ~ Z

    psdesign2 <- object + do.call(as.character(pdat[[1]]), pdat[-1])
    marg.est <- psdesign2 + ps_estimate()
    marg.VE <- mean(calc_risk(marg.est, contrast = "VE", bootstraps = FALSE)[, 2])
    emp.VE <- empirical_VE(object)
    cond.VE <- calc_risk(object, contrast = "VE", bootstraps = FALSE)
    cond.VE.est <- 1 - mean(cond.VE$R1)/mean(cond.VE$R0)
    VEtab <- c(empirical = emp.VE, marginal = marg.VE, model = cond.VE.est)

    empdiff <- 100 * VEtab[3]/VEtab[1] - 100
    mardiff <- 100 * VEtab[3]/VEtab[2] - 100

    cat("\nVaccine Efficacy: \n")
    print(VEtab, digits = digits)
    cat(sprintf("\tModel-based average VE is %.1f %% different from the empirical and %.1f %% different from the marginal.\n", empdiff, mardiff))
    if(abs(mardiff) > 100) warning("Check model and results carefully!")

    invisible(list(print = pout, VE.estimates = VEtab))

  }
}



#' Test for wide effect modificiation
#'
#' This runs a multivariate Wald test on the interaction terms of the model, using the bootstrap covariance
#'
#' @param x An object of class \link{psdesign} with bootstrap replicates
#'
#' @export
wem_test <- function(x){

  testdex <- grep(":Z", x$param.names, fixed = TRUE)
  bscov <- cov(x$bootstraps)

  if(length(testdex) == 1){

    bsse <- bscov[testdex, testdex]
    bsmu <- x$estimates$par[testdex]
    Xstat <- bsmu/bsse
    list(stat = Xstat, df = 1, p.value = 2 * pnorm(-abs(Xstat)))

  } else {

    Xstat <- x$estimates$par[testdex] %*% solve(bscov[testdex, testdex]) %*% x$estimates$par[testdex]
    list(stat = Xstat, df = length(testdex), p.value = pchisq(Xstat, df = length(testdex), lower.tail = FALSE))

  }

}