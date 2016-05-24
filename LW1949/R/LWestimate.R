#' Generate Litchfield and Wilcoxon Estimates
#'
#' Generate Litchfield and Wilcoxon's (1949) estimates in the evaluation of
#'   dose-effect experiments.
#' @param estParams
#'   A numeric vector of length 2, the estimated parameters (intercept and
#'   slope) of a straight line on the log10-probit scale, typically output
#'   from \code{\link{fitLWauto}}.
#' @param DEdata
#'   A data frame of dose-effect data (typically, the output from
#'   \code{\link{dataprep}}) containing at least eight variables:
#'   dose, ntot, nfx, pfx, log10dose, bitpfx, fxcateg, and LWkeep
#'   (see Details).
#' @details
#' The input data are expected to be summarized by dose.
#'   If duplicate doses are provided, an error will be thrown.
#' @return
#'   A list of length three:
#'   \itemize{
#'     \item \code{chi} = the chi-squared statistic with associated P value and
#'       degrees of freedom,
#'     \item \code{params} = the estimated intercept and slope of the
#'       dose-response curve on the log10 probit scale,
#'     \item \code{LWest} = the Litchfield Wilcoxon estimates of
#'       ED50 with 95\% confidence intervals and the number of records with
#'       partial effects (npartfx) as well as other metrics used in their
#'       step-by-step approach
#'       (ED16, ED84, S with 95\% confidence intervals, N', and fED50).
#'   }
#' @export
#' @import
#'   stats
#' @references
#'   Litchfield, JT Jr. and F Wilcoxon.  1949.
#'     A simplified method of evaluating dose-effect experiments.
#'     Journal of Pharmacology and Experimental Therapeutics 96(2):99-113.
#'     \href{http://jpet.aspetjournals.org/content/96/2/99.abstract}{[link]}.
#' @examples
#' dose <- c(0.0625, 0.125, 0.25, 0.5, 1)
#' ntested <- rep(8, 5)
#' nalive <- c(1, 4, 4, 7, 8)
#' mydat <- dataprep(dose=dose, ntot=ntested, nfx=nalive)
#' mydat
#' intslope <- fitLWauto(mydat)
#' LWestimate(intslope, mydat)

LWestimate <- function(estParams, DEdata) {
  if(length(estParams)!=2) stop("estParams must be of length 2")
  if(class(estParams)!="numeric") stop("estParams must be a numeric vector")
  if (!is.data.frame(DEdata)) stop("DEdata must be a data frame.")
  if (any(is.na(match(c("dose", "ntot", "nfx", "pfx", "log10dose", "bitpfx",
    "fxcateg", "LWkeep"), names(DEdata))))) {
    stop("DEdata must include these eight variables:",
      "dose, ntot, nfx, pfx, log10dose, bitpfx, fxcateg, LWkeep.")
  }
  dose.nona <- DEdata$dose[!is.na(DEdata$dose)]
  if (sum(duplicated(dose.nona))>0) {
    stop("Dose should be a vector of unique values, with no duplicates")
  }

  dfsub <- DEdata[DEdata$LWkeep, ]
  npartfx <- sum(dfsub$fxcateg==50)

  ### B1, B2, and C are all inside the function assessfit()
  # B1. If the expected value for any 0% or 100% dose is < 0.01% or > 99.99%,
  #   delete record
  # B2. Using the expected effect, record a corrected value for each
  #   0 and 100% effect
  ### C. The chi squared test
  chi <- assessfit(estParams, DEdata=dfsub, simple=FALSE)
  poorfit <- FALSE
  if (!is.na(chi$chi["pval"]) & chi$chi["pval"] < 0.05) {
    poorfit <- TRUE
    warning("Chi squared test indicates poor fit.")
  }

  ### D1. Read from the line on the graph the dose for 16, 50, and 84% effects
  ED16 <- predlinear(16, estParams, simple=TRUE)
  ED50 <- predlinear(50, estParams, simple=TRUE)
  ED84 <- predlinear(84, estParams, simple=TRUE)
  # D2. Calculate the Litchfield and Wilcoxon (1949) slope function, S
  S <- (ED84/ED50 + ED50/ED16) / 2
  # D3. Obtain Nprime, the total number of animals tested at those doses with
  #   expected effects between 16 and 84%.
  Nprime <- sum(dfsub$ntot[dfsub$dose >= ED16 & dfsub$dose <= ED84])

  # D4. Calculate S to the exponent for the ED50
  if(!poorfit) {
    fED50 <- S^(2.77/sqrt(Nprime))
  } else {
    # F. Factors for significantly heterogeneous data
    # value of t for p=0.05 and given df (Table 2)
    # n defined as df, K-2, in L-W's step C3
    n <- chi$chi["df"]
    t <- qt(1 - 0.05/2, n)
    formchunk <- as.vector(sqrt( chi$chi["chistat"] / (n*Nprime) ))
    fED50 <- S ^ ( 1.4 * t * formchunk )
  }

  # D5. Calculate the 95% confidence limits of the ED50
  if (is.finite(fED50)) {
    upper50 <- ED50 * fED50
    lower50 <- ED50 / fED50
  } else {
    upper50 <- NA
    lower50 <- NA
  }

  # E. Calculate the confidence limits of S
  #   (requires mysterious Nomograph No. 3)
  # E1. Calculate the dosage range as a ratio, R
  R <- max(dfsub$dose)/min(dfsub$dose)
  # E2. Calculate A from equation 6 in Appendix (Nomograph 3)
  A <- 10^( 1.1*(log10(S))^2 / log10(R) )

  # E3. Calculate K (the number of doses) and fS (Nomograph 2)
  K <- dim(dfsub)[1]
  if(!poorfit) {
    fS <- A^( 10*(K-1) / (K * sqrt(Nprime)) )
  } else {
    # F. Factors for significantly heterogeneous data
    fS = A ^ ( 5.1 * t * (K-1) * formchunk / K )
  }

  # E4. Calculate the 95% confidence limits of S
  upperS <- S * fS
  lowerS <- S / fS

  out <- list(chi=chi, params=estParams,
    LWest=c(ED50=ED50, lower=lower50, upper=upper50, npartfx=npartfx,
      ED16=ED16, ED84=ED84, S=S, lowerS=lowerS, upperS=upperS,
      Nprime=Nprime, fED50=fED50, fS=fS))
  return(out)
}
