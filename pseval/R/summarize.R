#' Calculate the risk and functions of the risk
#'
#' Computes the vaccince efficacy (VE) and other functions of the risk in each
#' treatment arm over the range of surrogate values observed in the data. VE(s)
#' is defined as 1 - risk(s, z = 1)/risk(s, z = 0), where z is the treatment
#' indicator. If any other variables are present in the risk model, then the
#' risk is computed at their median value.
#'
#' @details The contrast function is a function that takes 2 inputs, the risk_0
#' and risk_1, and returns some one dimensional function of those two inputs. It must be
#' vectorized. Some built-in functions are \code{"VE"} for vaccine efficacy = 1
#' - risk_1(s)/risk_0(s), \code{"RR"} for relative risk = risk_1(s)/risk_0(s),
#' \code{"logRR"} for log of the relative risk, and \code{"RD"} for the risk difference = risk_1(s) -
#' risk_0(s).
#'
#' @return A data frame containing columns for the S values, the computed contrast function at S, R0, and R1
#'   at those S values, and optionally standard errors and confidence intervals
#'   computed using bootstrapped estimates.
#'
#' @param psdesign A psdesign object. It must contain a risk model, an
#'   integration model, and estimated parameters. Bootstrapped parameters are
#'   optional
#' @param contrast The contrast function, or the name of the contrast function.
#'   See details.
#' @param t For time to event outcomes, a fixed time \code{t} may be provided to
#'   compute the cumulative distribution function. If not, the restricted mean
#'   survival time is used. Omit for binary outcomes.
#' @param sig.level Significance level for bootstrap confidence intervals
#' @param CI.type Character string, "pointwise" for pointwise confidence
#'   intervals, and "band" for simultaneous confidence band.
#' @param n.samps The number of samples to take over the range of S.1 at which
#'   the VE is calculated
#' @param bootstraps If true, and bootstrapped estimates are present, will
#'   calculate bootstrap standard errors and confidence bands.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # same result passing function name or function
#' calc_risk(binary.boot, contrast = "VE", n.samps = 20)
#' calc_risk(binary.boot, contrast = function(R0, R1) 1 - R1/R0, n.samps = 20)
#' }
calc_risk <- function(psdesign, contrast = "VE", t, sig.level = .05, CI.type = "band", n.samps = 5000, bootstraps = TRUE){

  stopifnot("estimates" %in% names(psdesign))

  impped <- psdesign$integration.models$S.1$icdf_sbarw(runif(n.samps))
  randrows <- sample(1:nrow(impped), ncol(impped), replace = TRUE)
  integrated <- impped[cbind(randrows, 1:n.samps)]
  obss <- psdesign$augdata$S.1

  trueobs <- sample(obss[!is.na(obss)],
                    floor(n.samps * mean(!is.na(obss))),
                    prob = psdesign$augdata$cdfweights[!is.na(obss)], replace = TRUE)


  if(is.factor(obss)){
    integrated <- factor(integrated, levels = levels(obss))
  }
  Splot <- sort(unlist(list(integrated, trueobs)))

  dat1 <- data.frame(S.1 = Splot, Z = 1)
  dat0 <- data.frame(S.1 = Splot, Z = 0)


  if(is.null(psdesign$risk.model$args$model)){

    others <- NULL

  } else {

    mod <- eval(psdesign$risk.model$args$model)
    mterms0 <- rownames(attr(terms(mod), "factors"))[-1]
    mterms <- unlist(lapply(colnames(psdesign$augdata), function(x){
      if(length(grep(x, mterms0, fixed = TRUE) > 0)){
        return(x)
      } else return(NULL)
    }))

    others <- mterms[!mterms %in% c("S.1", "Z")]
    for(j in others){

      tempval <- median(psdesign$augdata[, j], na.rm = TRUE)
      dat1[, j] <- tempval
      dat0[, j] <- tempval

    }

  }

  obsrisks <- riskcalc(psdesign$risk.function, psdesign$augdata$Y, psdesign$estimates$par,
                       t, dat0, dat1)

  obsVE <- data.frame(S.1 = Splot, Y = do.call(contrast, args = list(obsrisks$R0, obsrisks$R1)),
                      R0 = obsrisks$R0, R1 = obsrisks$R1)

  if(bootstraps && "bootstraps" %in% names(psdesign)){

    bsests <- psdesign$bootstraps
    bootYs <- matrix(NA, ncol = length(Splot) + 1, nrow = nrow(bsests))
    bootR0 <- matrix(NA, ncol = length(Splot) + 1, nrow = nrow(bsests))
    bootR1 <- matrix(NA, ncol = length(Splot) + 1, nrow = nrow(bsests))

    for(i in 1:nrow(bsests)){

      thispar <- as.numeric(bsests[i, -ncol(bsests)])
      thisrisks <- riskcalc(psdesign$risk.function, psdesign$augdata$Y, thispar,
                       t, dat0, dat1)

      bootYs[i, ] <- c(do.call(contrast, list(thisrisks$R0, thisrisks$R1)), bsests[i, "convergence"])
      bootR0[i, ] <- c(thisrisks$R0, bsests[i, "convergence"])
      bootR1[i, ] <- c(thisrisks$R1, bsests[i, "convergence"])

    }

    bootYs <- as.data.frame(bootYs)
    bootR1 <- as.data.frame(bootR1)
    bootR0 <- as.data.frame(bootR0)

    colnames(bootYs)[ncol(bootYs)] <- colnames(bootR0)[ncol(bootR0)] <- colnames(bootR1)[ncol(bootR1)] <- "convergence"
    A1 <- as.data.frame(summarize_bs(bootR1, obsVE$R1, sig.level = sig.level, CI.type = CI.type)$table)
    A2 <- as.data.frame(summarize_bs(bootR0, obsVE$R0, sig.level = sig.level, CI.type = CI.type)$table)
    A3 <- as.data.frame(summarize_bs(bootYs, obsVE$Y, sig.level = sig.level, CI.type = CI.type)$table)

    colnames(A1) <- gsub("%", "", paste("R1", colnames(A1), sep = "."), fixed = TRUE)
    colnames(A2) <- gsub("%", "", paste("R0", colnames(A2), sep = "."), fixed = TRUE)
    colnames(A3) <- gsub("%", "", paste("Y", colnames(A3), sep = "."), fixed = TRUE)

    obsVE <- cbind(obsVE, A3, A2, A1)

  }

  attr(obsVE, "Y.function") <- substitute(contrast)
  obsVE

}



#' Summarize bootstrap samples
#'
#' @param bootdf Data frame containing bootstrapped estimates, with a column containing a convergence indicator
#' @param estdf Data frame containing full sample estimate
#' @param sig.level Significance level to use for confidence intervals
#' @param CI.type Character string, "pointwise" for pointwise confidence intervals, and "band" for simultaneous confidence band.
#' @export

summarize_bs <- function(bootdf, estdf = NULL, sig.level = .05, CI.type = "band") {

  bs <- bootdf[bootdf$convergence == 0, -which(colnames(bootdf) == "convergence"), drop = FALSE]

  if(CI.type == "pointwise"){
    mary <- function(x){
      funlist <- list(boot.se = stats::sd,
                    lower.CL = function(x) quantile(x, sig.level/2),
                    upper.CL = function(x) quantile(x, 1 - sig.level/2))

      sapply(funlist, function(f) f(x))
    }

    table <- t(sapply(bs, mary))

  } else if(CI.type == "band"){

    estmat <- matrix(rep(estdf, nrow(bs)), nrow = nrow(bs), byrow = TRUE)
    maxdiff <- apply(abs(as.matrix(bs) - estmat), MARGIN = 1, FUN = max)
    inCI <- as.matrix(bs)[which(maxdiff < quantile(maxdiff, 1 - sig.level)), ]
    upper.CL <- apply(inCI, MARGIN = 2, FUN = max)
    lower.CL <- apply(inCI, MARGIN = 2, FUN = min)

    table0 <- data.frame(boot.se = sapply(bs, stats::sd))
    table <- data.frame(upper.CL = upper.CL, lower.CL = lower.CL)
    colnames(table) <- paste(colnames(table), 1-sig.level, sep = ".")

    table <- cbind(table0, table)

  }

  conv <- c(nboot = nrow(bootdf), ncov = sum(bootdf$convergence == 0))

  list(table = table, conv = conv)

}


#' Compute the empirical Vaccine Efficacy
#'
#' @param psdesign An object of class \link{psdesign}
#' @param t Fixed time for time to event outcomes to compute VE. If missing, uses restricted mean survival.
#'
#' @export
empirical_VE <- function(psdesign, t){

  pd <- psdesign$augdata
  if(inherits(pd$Y, "Surv")){

    if(missing(t)){

      ttt <- summary(survival::survfit(pd$Y ~ 1), rmean = "common")$table[["*rmean"]]

    } else ttt <- t

    sf <- survival::survfit(pd$Y ~ pd$Z)
    sfs <- 1 - summary(sf, times = ttt, extend = TRUE)$surv

    1 - sfs[2]/sfs[1]

  } else {

    1 - mean(pd$Y[pd$Z == 1])/mean(pd$Y[pd$Z == 0])

  }

}


#' Vaccine efficacy contrast functions
#'
#' @param R0 A vector of risks in the control arm
#' @param R1 A vector of risks in the treatment arm
#'
#' @keywords Internal
#'
#' @return A vector the same length as R0 and R1.
#'
#' @details These functions take the risk in the two treatment arms, and
#'   computes a one-dimensional summary of those risks. Built-in choices are
#'   \code{"VE"} for vaccine efficacy = 1 - risk_1(s)/risk_0(s), \code{"RR"} for
#'   relative risk = risk_0(s)/risk_1(s), \code{"logRR"} for log of the relative
#'   risk, and \code{"RD"} for the risk difference = risk_0(s) - risk_1(s).

VE <- function(R0, R1){
  1 - R1/R0
}

RR <- function(R0, R1){
  R1/R0
}

logRR <- function(R0, R1){
  log(R1/R0)
}

RD <- function(R0, R1){
  R1 - R0
}

#' Calculate risks with handlers for survival data
#' @keywords Internal
#' @param risk.function Function taking three arguments, a data.frame, parameters, and time. It should return a vector the same number of rows as the data frame
#' @param Y The outcome variable
#' @param par the vector of parameter values
#' @param t Time for a survival outcome, may be missing
#' @param dat0 Data frame containing S and Z = 1
#' @param dat1 Data frame containing S and Z = 0
#'

riskcalc <- function(risk.function, Y, par, t, dat0, dat1){

    if(inherits(Y, "Surv") && missing(t)){

    ttt <- summary(survival::survfit(Y ~ 1), rmean = "common")$table[["*rmean"]]

    warning(sprintf("No time given for time to event outcome, using restricted mean survival: %.1f", ttt))
    R1 <- risk.function(dat1, par, t = ttt)
    R0 <- risk.function(dat0, par, t = ttt)


  } else if(inherits(Y, "Surv")) {

    R1 <- risk.function(dat1, par, t)
    R0 <- risk.function(dat0, par, t)

  } else {

    R1 <- risk.function(dat1, par)
    R0 <- risk.function(dat0, par)

  }

  list(R0 = R0, R1 = R1)

}