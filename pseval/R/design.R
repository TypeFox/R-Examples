#' Specify a design for a principal surrogate evaluation
#'
#' Generate mappings that describe how variables in the data are mapped to
#' components of the principal surrogate analysis. Other than \code{data}, this
#' is a list of key-value pairs desribing the common elements of a ps analysis.
#' The required keys are Z, Y, and S. Optional keys are BIP, CPV, BSM, and
#' weights. These elements are described in details below. Additional keys-value
#' pairs can be included in \code{...}. This function generates an augmented
#' dataset and additional information for subsequent steps in the analysis. In
#' the subsequent steps, refer to the variables by the keys. See
#' \link{add_integration} and \link{add_riskmodel} for information on how to
#' proceed in the analysis.
#'
#' @param data Data frame containing data to be analyzed
#' @param Z Expression defining the treatment variable which has 2 levels
#' @param Y Expression defining the outcome variable. For binary events this
#'   should be coded as 0/1 or a factor with 2 levels. For censored
#'   time-to-event outcomes this can be a call to \link[survival]{Surv}
#' @param S Expression defining the candidate surrogate
#' @param BIP Optional expression defining the baseline irrelevant predictor
#' @param CPV Optional expression defining the closeout placebo vaccination
#'   measurement
#' @param BSM Optional expression defining the baseline surrogate measurement
#' @param weights optional expression defining weights to accomodate nonrandom
#'   subsampling, such as case control or two phase
#' @param tau numeric, When the outcome Y is a survival time, it is possible
#'   that the surrogate was measured at some time tau after enrollment. Use the
#'   argument tau to specify the time when the surrogate was measured, in study
#'   time. Not required for binary Y.
#' @param ... Other key-value pairs that will be included in the augmented data,
#'   e.g. additional candidate surrogates, covariates for adjustment, variables
#'   used for integration
#'
#' @export
#' @importFrom graphics lines plot segments
#' @importFrom stats cov dpois gaussian glm median model.frame model.matrix model.offset model.response optim pchisq pnorm ppois predict qnorm quantile quasi quasibinomial rbinom rexp rnorm runif terms
#' @importFrom utils flush.console head setTxtProgressBar txtProgressBar

psdesign <- function(data, Z, Y, S,
                     BIP = NULL, CPV = NULL, BSM = NULL, weights = NULL, tau, ...){


  mapping <- sapply(as.list(match.call())[-c(1, 2)], as.character)

  ## if Y is a survival object, check for presence of tau, the time when the biomarker was measured

  oot <- eval(substitute(Y), envir = data)
  if(inherits(oot, "Surv")){

    if(missing(tau)){
      warning("tau missing in psdesign: assuming that the surrogate S was measured at time 0.")
    } else {

      oot[, 1] <- oot[, 1] - tau
      nminus <- sum(oot[, 1] < 0)

      stdex <- oot[, 1] >= 0
      if(nminus > 0){
        warning("Removing ", nminus, " subjects with survival times less that tau.")
      }
      data <- data[stdex, ]

    }

  }

  trt <- verify_trt(eval(substitute(Z), envir = data))

  oot <- eval(substitute(Y), envir = data)
  #if(!inherits(oot, "Surv")){
  #  oot <- verify_trt(oot)
  #}

  cand <- eval(substitute(S), envir = data)

  S.0 <- S.1 <- cand
  S.0[trt != 0] <- NA
  S.1[trt != 1] <- NA
  ## some logic based on presence of BIP/CPV/BSM

  bip <- eval(substitute(BIP), envir = data)
  cpv <- eval(substitute(CPV), envir = data)
  bsm <- eval(substitute(BSM), envir = data)

  if(!is.null(cpv)){
    S.1[is.na(S.1) & no_event(oot)] <- cpv[is.na(S.1) & no_event(oot)]
  }
  if(!is.null(bsm)){
    S.0[is.na(S.0)] <- bsm[is.na(S.0)]
  }

  if(is.null(weights)){
    weights <- rep(1, length(S.0))
  } else weights <- eval(substitute(weights), envir = data)
  rval <- NULL

  rval$augdata <- data.frame(Z = trt, Y = oot, S.1 = S.1, S.0 = S.0,
                             cdfweights = weights)

  if(!is.null(bip)) rval$augdata <- data.frame(rval$augdata, BIP = bip)
  if(!is.null(cpv)) rval$augdata <- data.frame(rval$augdata, CPV = cpv)
  if(!is.null(bsm)) rval$augdata <- data.frame(rval$augdata, BSM = bsm)


  optionals <- do.call(cbind, eval(substitute(list(...)), envir = data))
  if(!is.null(optionals)) rval$augdata <- data.frame(rval$augdata, optionals)

  rval$mapping <- mapping

  class(rval) <- c("ps", "psdesign")

  rval

}


#' Check that a variable is suitable for using as binary treatment indicator
#'
#' Checks for two classes and gives a warning message indicating which level is assumed to be 0/1
#'
#' @param D Vector that will be checked for 2-class labels
#'
#' @keywords Internal
#'
verify_trt <- function(D){

  if(length(levels(as.factor(D))) > 2) stop("Only variables with 2 levels supported")

  slev <- sort(levels(as.factor(D)))
  if(slev[1] == 0 & slev[2] == 1) return(D)

  warning(paste0("Variable not labeled 0/1, assuming ", slev[1], " = 0 and ", slev[2], " = 1!"))

  zero1 <- c(0, 1)
  names(zero1) <- slev

  zero1[D]

}


no_event <- function(Y){

  if(inherits(Y, "Surv")){
    stopifnot(dim(Y)[2] <= 2)
    Y[, 2] == 0
  } else {
    Y == 0
  }
}

