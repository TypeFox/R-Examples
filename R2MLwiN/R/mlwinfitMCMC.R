#' An S4 class that stores the outputs of the fitted MCMC model.
#'
#' An MLwiN model run via the MCMC estimation method is represented by an "mlwinfitMCMC" object
#'
#' @section An instance of the Class:
#'  An instance is created by calling function \code{\link{runMLwiN}}.
#'
#' @slot Nobs Computes the number of complete observations.
#' @slot DataLength Total number of cases.
#' @slot Hierarchy For each higher level of a multilevel model, returns the number of units at that level, together with the minimum, mean and maximum number of lower-level units nested within units of the current level.
#' @slot burnin An integer specifying length of the burn-in.
#' @slot nchains An integer specifying number of MCMC chains run.
#' @slot iterations An integer specifying the number of iterations after burn-in.
#' @slot D A vector specifying the type of distribution to be modelled, which can include \code{'Normal'}, \code{'Binomial'} \code{'Poisson'}, \code{'Multinomial'}, \code{'Multivariate Normal'}, or \code{'Mixed'}.
#' @slot Formula A formula object (or a character string) specifying a multilevel model.
#' @slot levID A character string (vector) of the specified level ID(s).
#' @slot merr A vector which sets-up measurement errors on predictor variables.
#' @slot fact A list of objects specified for factor analysis, including \code{nfact}, \code{lev.fact}, \code{nfactor}, \code{factor}, \code{loading} and \code{constr}.
#' @slot xc A list of objects specified for cross-classified and/or multiple membership models, including \code{class}, \code{N1}, \code{weight}, \code{id} and \code{car}. 
#' @slot FP Displays the fixed part estimates.
#' @slot RP Displays the random part estimates.
#' @slot FP.cov Displays a covariance matrix of the fixed part estimates.
#' @slot RP.cov Displays a covariance matrix of the random part estimates.
#' @slot chains Captures the MCMC chains from MLwiN for all parameters.
#' @slot elapsed.time Calculates the CPU time used for fitting the model.
#' @slot BDIC Bayesian Deviance Information Criterion (DIC)
#' @slot call The matched call.
#' @slot LIKE The deviance statistic (-2*log(like)).
#' @slot fact.loadings If \code{fact} is not empty, then the factor loadings are returned.
#' @slot fact.loadings.sd If \code{fact} is not empty, then the factor loading standard deviationss are returned.
#' @slot fact.cov If \code{fact} is not empty, then factor covariances are returned.
#' @slot fact.cov.sd If \code{fact} is not empty, then factor covariance standard deviations are returned.
#' @slot fact.chains If \code{fact} is not empty, then the factor chains are returned.
#' @slot MIdata If \code{dami[1]} is one then the mean complete response variable \code{y} is returned for each chain, if \code{dami[1]} is two then the SD is also included.
#' @slot imputations If \code{dami[1]} is zero, then a list of completed datasets containing complete response variable \code{y} is returned.
#' @slot residual If \code{resi.store} is \code{TRUE}, then the residual estimates at all levels are returned.
#' @slot resi.chains If \code{resi.store.levs} is not empty, then the residual chains at these levels are returned.
#' @slot version The MLwiN version used to fit the model
#' @slot data The data.frame that was used to fit the model.
#'
#' @author Zhang, Z., Charlton, C.M.J., Parker, R.M.A., Leckie, G., and Browne,
#' W.J. (2015) Centre for Multilevel Modelling, University of Bristol.
#'
#' @seealso
#' \code{\link{runMLwiN}}
#'
#' @examples
#' \dontrun{
#' 
#' library(R2MLwiN)
#' # NOTE: if MLwiN not saved in location R2MLwiN defaults to, specify path via:
#' # options(MLwiN_path = 'path/to/MLwiN vX.XX/')
#' # If using R2MLwiN via WINE, the path may look like this:
#' # options(MLwiN_path = '/home/USERNAME/.wine/drive_c/Program Files (x86)/MLwiN vX.XX/')
#'   
#' ## Example: tutorial
#' data(tutorial, package = "R2MLwiN")
#'
#' (mymodel <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school) + (1 | student),
#'                      estoptions = list(EstM = 1), data = tutorial)) 
#' 
#' ##summary method
#' summary(mymodel)
#' 
#' ##BDIC slot
#' mymodel@@BDIC
#' }
#'
#' @name mlwinfitMCMC-class
#' @rdname mlwinfitMCMC-class
#' @exportClass mlwinfitMCMC
setClass(Class = "mlwinfitMCMC", representation = representation(version = "character", Nobs = "numeric", DataLength = "numeric", Hierarchy = "ANY",
                                                                 burnin = "numeric", iterations = "numeric", nchains = "numeric", D = "ANY", Formula = "ANY", levID = "character", 
                                                                 merr = "ANY", fact = "ANY", xc = "ANY", FP = "numeric", RP = "numeric", RP.cov = "matrix", FP.cov = "matrix", 
                                                                 chains = "ANY", elapsed.time = "numeric", call = "ANY", BDIC = "numeric", LIKE = "ANY", fact.loadings = "numeric", 
                                                                 fact.loadings.sd = "numeric", fact.cov = "numeric", fact.cov.sd = "numeric", fact.chains = "ANY", MIdata = "data.frame", 
                                                                 imputations = "list", residual = "list", resi.chains = "ANY", data = "data.frame"))

#' Extract or Replace parts of "mlwinfitMCMC" objects
#' @param x data frame
#' @param i,j elements to extract or replace. For \code{[} and \code{[[}, these are \code{character}.
#' @param drop not used.
#' @param value a suitable replacement value.
#' @rdname extract-methods-mcmc
setMethod("[", "mlwinfitMCMC", function(x, i, j, drop) {
  if (i == "version") {
    return(x@version)
  }
  if (i == "Nobs") {
    return(x@Nobs)
  }
  if (i == "DataLength") {
    return(x@DataLength)
  }
  if (i == "Hierarchy") {
    return(x@Hierarchy)
  }
  if (i == "burnin") {
    return(x@burnin)
  }
  if (i == "iterations") {
    return(x@iterations)
  }
  if (i == "nchains") {
    return(x@nchains)
  }
  if (i == "D") {
    return(x@D)
  }
  if (i == "Formula") {
    return(x@Formula)
  }
  if (i == "levID") {
    return(x@levID)
  }
  if (i == "merr") {
    return(x@merr)
  }
  if (i == "fact") {
    return(x@fact)
  }
  if (i == "xc") {
    return(x@xc)
  }
  if (i == "FP") {
    return(x@FP)
  }
  if (i == "RP") {
    return(x@RP)
  }
  if (i == "FP.cov") {
    return(x@FP.cov)
  }
  if (i == "RP.cov") {
    return(x@RP.cov)
  }
  if (i == "chains") {
    return(x@chains)
  }
  if (i == "elapsed.time") {
    return(x@elapsed.time)
  }
  if (i == "BDIC") {
    return(x@BDIC)
  }
  if (i == "call") {
    return(x@call)
  }
  if (i == "LIKE") {
    return(x@LIKE)
  }
  if (i == "fact.loadings") {
    return(x@fact.loadings)
  }
  if (i == "fact.loadings.sd") {
    return(x@fact.loadings.sd)
  }
  if (i == "fact.cov") {
    return(x@fact.cov)
  }
  if (i == "fact.cov.sd") {
    return(x@fact.cov.sd)
  }
  if (i == "fact.chains") {
    return(x@fact.chains)
  }
  if (i == "MIdata") {
    return(x@MIdata)
  }
  if (i == "imputations") {
    return(x@imputations)
  }
  if (i == "residual") {
    return(x@residual)
  }
  if (i == "resi.chains") {
    return(x@resi.chains)
  }
  if (i == "data") {
    return(x@data)
  }
})

#' @rdname extract-methods-mcmc
setReplaceMethod("[", signature(x = "mlwinfitMCMC"), function(x, i, j, value) {
  if (i == "version") {
    x@version <- value
  }
  if (i == "Nobs") {
    x@Nobs <- value
  }
  if (i == "DataLength") {
    x@DataLength <- value
  }
  if (i == "Hierarchy") {
    x@Hierarchy <- value
  }
  if (i == "burnin") {
    x@burnin <- value
  }
  if (i == "iterations") {
    x@iterations <- value
  }
  if (i == "nchains") {
    x@nchains <- value
  }
  if (i == "D") {
    x@D <- value
  }
  if (i == "Formula") {
    x@Formula <- value
  }
  if (i == "levID") {
    x@levID <- value
  }
  if (i == "merr") {
    x@merr <- value
  }
  if (i == "fact") {
    x@fact <- value
  }
  if (i == "xc") {
    x@xc <- value
  }
  if (i == "FP") {
    x@FP <- value
  }
  if (i == "RP") {
    x@RP <- value
  }
  if (i == "FP.cov") {
    x@FP.cov <- value
  }
  if (i == "RP.cov") {
    x@RP.cov <- value
  }
  if (i == "chains") {
    x@chains <- value
  }
  if (i == "elapsed.time") {
    x@elapsed.time <- value
  }
  if (i == "BDIC") {
    x@BDIC <- value
  }
  if (i == "call") {
    x@call <- value
  }
  if (i == "LIKE") {
    x@LIKE <- value
  }
  if (i == "fact.loadings") {
    x@fact.loadings <- value
  }
  if (i == "fact.loadings.sd") {
    x@fact.loadings.sd <- value
  }
  if (i == "fact.cov") {
    x@fact.cov <- value
  }
  if (i == "fact.cov.sd") {
    x@fact.cov.sd <- value
  }
  if (i == "fact.chains") {
    x@fact.chains <- value
  }
  if (i == "MIdata") {
    x@MIdata <- value
  }
  if (i == "imputations") {
    x@imputations <- value
  }
  if (i == "residual") {
    x@residual <- value
  }
  if (i == "resi.chains") {
    x@resi.chains <- value
  }
  if (i == "data") {
    x@data <- value
  }
  validObject(x)
  return(x)
})

#' @rdname extract-methods-mcmc
setMethod("[[", "mlwinfitMCMC", function(x, i, j, drop) {
  if (i == "version") {
    return(x@version)
  }
  if (i == "Nobs") {
    return(x@Nobs)
  }
  if (i == "DataLength") {
    return(x@DataLength)
  }
  if (i == "Hierarchy") {
    return(x@Hierarchy)
  }
  if (i == "burnin") {
    return(x@burnin)
  }
  if (i == "iterations") {
    return(x@iterations)
  }
  if (i == "nchains") {
    return(x@nchains)
  }
  if (i == "D") {
    return(x@D)
  }
  if (i == "Formula") {
    return(x@Formula)
  }
  if (i == "levID") {
    return(x@levID)
  }
  if (i == "merr") {
    return(x@merr)
  }
  if (i == "fact") {
    return(x@fact)
  }
  if (i == "xc") {
    return(x@xc)
  }
  if (i == "FP") {
    return(x@FP)
  }
  if (i == "RP") {
    return(x@RP)
  }
  if (i == "FP.cov") {
    return(x@FP.cov)
  }
  if (i == "RP.cov") {
    return(x@RP.cov)
  }
  if (i == "chains") {
    return(x@chains)
  }
  if (i == "elapsed.time") {
    return(x@elapsed.time)
  }
  if (i == "BDIC") {
    return(x@BDIC)
  }
  if (i == "call") {
    return(x@call)
  }
  if (i == "LIKE") {
    return(x@LIKE)
  }
  if (i == "fact.loadings") {
    return(x@fact.loadings)
  }
  if (i == "fact.loadings.sd") {
    return(x@fact.loadings.sd)
  }
  if (i == "fact.cov") {
    return(x@fact.cov)
  }
  if (i == "fact.cov.sd") {
    return(x@fact.cov.sd)
  }
  if (i == "fact.chains") {
    return(x@fact.chains)
  }
  if (i == "MIdata") {
    return(x@MIdata)
  }
  if (i == "imputations") {
    return(x@imputations)
  }
  if (i == "residual") {
    return(x@residual)
  }
  if (i == "resi.chains") {
    return(x@resi.chains)
  }
  if (i == "data") {
    return(x@data)
  }
})

#' @rdname extract-methods-mcmc
setReplaceMethod("[[", signature(x = "mlwinfitMCMC"), function(x, i, j, value) {
  if (i == "version") {
    x@version <- value
  }
  if (i == "Nobs") {
    x@Nobs <- value
  }
  if (i == "DataLength") {
    x@DataLength <- value
  }
  if (i == "Hierarchy") {
    x@Hierarchy <- value
  }
  if (i == "burnin") {
    x@burnin <- value
  }
  if (i == "iterations") {
    x@iterations <- value
  }
  if (i == "nchains") {
    x@nchains <- value
  }
  if (i == "D") {
    x@D <- value
  }
  if (i == "Formula") {
    x@Formula <- value
  }
  if (i == "levID") {
    x@levID <- value
  }
  if (i == "merr") {
    x@merr <- value
  }
  if (i == "fact") {
    x@fact <- value
  }
  if (i == "xc") {
    x@xc <- value
  }
  if (i == "FP") {
    x@FP <- value
  }
  if (i == "RP") {
    x@RP <- value
  }
  if (i == "FP.cov") {
    x@FP.cov <- value
  }
  if (i == "RP.cov") {
    x@RP.cov <- value
  }
  if (i == "chains") {
    x@chains <- value
  }
  if (i == "elapsed.time") {
    x@elapsed.time <- value
  }
  if (i == "BDIC") {
    x@BDIC <- value
  }
  if (i == "call") {
    x@call <- value
  }
  if (i == "LIKE") {
    x@LIKE <- value
  }
  if (i == "fact.loadings") {
    x@fact.loadings <- value
  }
  if (i == "fact.loadings.sd") {
    x@fact.loadings.sd <- value
  }
  if (i == "fact.cov") {
    x@fact.cov <- value
  }
  if (i == "fact.cov.sd") {
    x@fact.cov.sd <- value
  }
  if (i == "fact.chains") {
    x@fact.chains <- value
  }
  if (i == "MIdata") {
    x@MIdata <- value
  }
  if (i == "imputations") {
    x@imputations <- value
  }
  if (i == "residual") {
    x@residual <- value
  }
  if (i == "resi.chains") {
    x@resi.chains <- value
  }
  if (i == "data") {
    x@data <- value
  }
  validObject(x)
  return(x)
})

#' Summarize "mlwinfitMCMC" objects
#' @param object,x an \code{\link{mlwinfitMCMC-class}} object
#' @param ... other parameters
#' @param digits the number of significant digits to use when printing.
#' @param signif.stars logical. If TRUE, 'significance stars' are printed for each coefficient.
#' @param z.ratio logical. If TRUE, z-ratio values are displayed for each coefficient.
#' @rdname summary-methods-mcmc
#' @export 
setMethod("summary", signature(object = "mlwinfitMCMC"), function(object, ...) {
  object
})

printMCMC <- function(x, digits = max(3, getOption("digits") - 2), signif.stars = getOption("show.signif.stars"), 
                      z.ratio = TRUE, ...) {
  
  object <- summary(x)
  align2right <- function(titlename, ele) {
    # for printing the table on the screen
    all.ele <- c(titlename, ele)
    len.all.ele <- nchar(all.ele)
    max.len.ele <- max(len.all.ele)
    for (j in 1:length(all.ele)) {
      if (len.all.ele[j] < max.len.ele) {
        len.diff <- max.len.ele - len.all.ele[j]
        all.ele[j] <- paste(paste(rep(" ", len.diff), collapse = ""), all.ele[j], sep = "")
      }
    }
    all.ele
  }
  
  align2left <- function(titlename, ele) {
    # for printing the table on the screen
    all.ele <- c(titlename, ele)
    len.all.ele <- nchar(all.ele)
    max.len.ele <- max(len.all.ele)
    for (j in 1:length(all.ele)) {
      if (len.all.ele[j] < max.len.ele) {
        len.diff <- max.len.ele - len.all.ele[j]
        all.ele[j] <- paste(all.ele[j], paste(rep(" ", len.diff), collapse = ""), sep = "")
      }
    }
    all.ele
  }
  
  signifstar <- function(pval) {
    starstr <- "N/A"
    if (!is.na(pval) && pval >= 0 && pval <= 1) {
      if (pval < 0.001) {
        starstr <- "***"
      }
      if (pval >= 0.001 && pval < 0.01) {
        starstr <- "** "
      }
      if (pval >= 0.01 && pval < 0.05) {
        starstr <- "*  "
      }
      if (pval >= 0.05 && pval < 0.1) {
        starstr <- ".  "
      }
      if (pval >= 0.1) {
        starstr <- "   "
      }
    }
    starstr
  }
  
  chainnames <- coda::varnames(object@chains)
  FP.names <- grep("^FP_", chainnames, value = TRUE)
  RP.names <- grep("^RP[0-9]+_", chainnames, value = TRUE)
  ESS <- effectiveSize(object@chains)
  levID0 <- object@levID
  cat("\n")
  cat(paste(rep("-", 50), collapse = "*"), "\n")
  cat(object@version, " multilevel model", paste("(", object@D[1], ")", sep = ""), "\n")
  if (!is.null(object@Hierarchy)) 
    print(object@Hierarchy)
  if (!is.null(object@xc) && isTRUE(object@xc)) {
    cat("Estimation algorithm:  MCMC      Cross-classified              Elapsed time :", paste(round(object@elapsed.time, 
                                                                                                     2), "s", sep = ""), "\n")
  } else {
    cat("Estimation algorithm:  MCMC      Elapsed time :", paste(round(object@elapsed.time, 2), "s", sep = ""), 
        "\n")
  }
  cat("Number of obs: ", object@Nobs, paste0("(from total ", object@DataLength, ")"), "         Number of iter.:", object@iterations, 
      " Chains:", object@nchains, " Burn-in:", object@burnin, "\n")
  
  if (!(object@D[1] == "Mixed") && is.null(object@merr) && is.null(object@fact)) {
    cat("Bayesian Deviance Information Criterion (DIC)\n")
    cat("Dbar      D(thetabar)    pD      DIC\n")
    cat(formatC(object@BDIC, format = "f", digits = 3, width = -10), "\n")
  } else {
    cat(paste("Deviance statistic: ", round(object@LIKE, 1)), "\n")
  }
  
  
  cat(paste(rep("-", 50), collapse = "-"), "\n")
  cat("The model formula:\n")
  print(object@Formula)
  levID.display <- ""
  if (is.na(levID0[length(levID0)])) {
    levID0 <- levID0[-length(levID0)]
  }
  for (i in 1:length(levID0)) {
    levID.display <- paste(levID.display, "Level ", length(levID0) + 1 - i, ": ", levID0[i], "     ", sep = "")
  }
  cat(levID.display, "\n")
  cat(paste(rep("-", 50), collapse = "-"), "\n")
  
  if (!is.null(object@fact) && object@D[1] == "Multivariate Normal") {
    qt025 <- object@fact.loadings - qnorm(0.975) * object@fact.loadings.sd
    qt975 <- object@fact.loadings + qnorm(0.975) * object@fact.loadings.sd
    loads <- rbind(object@fact.loadings, object@fact.loadings.sd, qt025, qt975)
    
    for (j in 1:object@fact$nfact) {
      cat("The estimates of factor", j, "loadings:\n")
      loadx.names <- colnames(loads)[grep(paste0("load+", j, "+\\_"), colnames(loads))]
      loadx <- loads[, loadx.names]
      printcol0 <- align2left("        ", loadx.names)
      printcol1 <- align2right("Coef.", format(round(loadx[1, ], digits), nsmall = digits))
      printcol2 <- align2right("Std. Err.", format(round(loadx[2, ], digits), nsmall = digits))
      printcol3 <- align2right("[95% Conf.", format(round(loadx[3, ], digits), nsmall = digits))
      printcol4 <- align2right("Interval]", format(round(loadx[4, ], digits), nsmall = digits))
      for (i in 1:(ncol(loadx) + 1)) {
        cat(printcol0[i], " ", printcol1[i], " ", printcol2[i], " ", printcol3[i], " ", printcol4[i], "\n")
      }
      cat(paste(rep("-", 50), collapse = "-"), "\n")
    }
    
    qt025 <- object@fact.cov - qnorm(0.975) * object@fact.cov.sd
    qt975 <- object@fact.cov + qnorm(0.975) * object@fact.cov.sd
    fcov <- rbind(object@fact.cov, object@fact.cov.sd, qt025, qt975)
    
    cat("The estimates of factor covariances:\n")
    printcol0 <- align2left("        ", colnames(fcov))
    printcol1 <- align2right("Coef.", format(round(fcov[1, ], digits), nsmall = digits))
    printcol2 <- align2right("Std. Err.", format(round(fcov[2, ], digits), nsmall = digits))
    printcol3 <- align2right("[95% Conf.", format(round(fcov[3, ], digits), nsmall = digits))
    printcol4 <- align2right("Interval]", format(round(fcov[4, ], digits), nsmall = digits))
    for (i in 1:(ncol(fcov) + 1)) {
      cat(printcol0[i], " ", printcol1[i], " ", printcol2[i], " ", printcol3[i], " ", printcol4[i], "\n")
    }
    cat(paste(rep("-", 50), collapse = "-"), "\n")
  }
  
  cat("The fixed part estimates: ", "\n")

  chain.stats <- summary(object@chains, quantiles=c(0.025, 0.975))
  chain.means <- chain.stats$statistics[,"Mean"]
  chain.sds <- chain.stats$statistics[,"SD"]
  chain.qt025 <- chain.stats$quantiles[,"2.5%"]
  chain.qt975 <- chain.stats$quantiles[,"97.5%"]

  if (sum(grepl("bcons", colnames(object@chains))) > 0) {
    bcons.pos <- grep("bcons", colnames(object@chains))
    object@chains[1, bcons.pos] <- object@chains[1, bcons.pos] - 0.001
  }

  t.stats <- chain.means / chain.sds
  
  p.values <- 2 * pnorm(abs(t.stats), lower.tail = FALSE)
  t.stat <- NULL
  for (i in FP.names) t.stat <- c(t.stat, t.stats[[i]])
  p.value <- NULL
  for (i in FP.names) p.value <- c(p.value, p.values[[i]])
  onesided.p.value <- NULL
  for (i in FP.names) {
    x <- as.matrix(object@chains[, i])
    if (sign(mean(x)) > 0) {
      onesided.p.values <- sum(x < 0)/length(x)
    } else {
      onesided.p.values <- sum(x > 0)/length(x)
    }
    onesided.p.value <- c(onesided.p.value, onesided.p.values)
  }
  strstar <- as.vector(sapply(p.value, signifstar))
  FP.print <- rbind(chain.means[FP.names], chain.sds[FP.names], t.stat, p.value, onesided.p.value, chain.qt025[FP.names], chain.qt975[FP.names], ESS[FP.names])
  FP.names2 <- gsub("FP+\\_", "", FP.names)
  
  printcol0 <- align2left("        ", FP.names2)
  printcol1 <- align2right("Coef.", format(round(FP.print[1, ], digits), nsmall = digits))
  printcol2 <- align2right("Std. Err.", format(round(FP.print[2, ], digits), nsmall = digits))
  printcol3 <- align2right("z", format(round(FP.print[3, ], 2), nsmall = 2))
  printcol4 <- align2right("Pr(>|z|)", formatC(FP.print[4, ]))
  printcol4b <- align2right("   ", strstar)
  printcol5 <- align2right("pMCMC(1-sided)", formatC(FP.print[5, ]))
  printcol6 <- align2right("[95% Cred.", format(round(FP.print[6, ], digits), nsmall = digits))
  printcol7 <- align2right("Interval]", format(round(FP.print[7, ], digits), nsmall = digits))
  printcol8 <- align2right("ESS", format(round(FP.print[8, ]), nsmall = 0))
  for (i in 1:(ncol(FP.print) + 1)) {
    cat(printcol0[i], " ", printcol1[i], " ", printcol2[i], " ")
    if (z.ratio) {
      cat(printcol3[i], " ", printcol4[i], " ")
      if (signif.stars) {
        cat(printcol4b[i], " ")
      }
    } else {
      cat(printcol5[i], " ")
    }
    cat(printcol6[i], " ", printcol7[i], " ", printcol8[i], "\n")
  }
  if (signif.stars & z.ratio) {
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ", "\n")
  }
  
  nlev <- length(object@levID)
  if (is.na(object@levID[length(object@levID)])) {
    mlwinlev <- (nlev - 1):1
    levID2 <- levID0
  } else {
    mlwinlev <- nlev:1
    levID2 <- object@levID
  }
  
  RP.print <- rbind(chain.means[RP.names], chain.sds[RP.names], chain.qt025[RP.names], chain.qt975[RP.names], ESS[RP.names])
  for (i in 1:length(mlwinlev)) {
    RPx.pos <- grep(paste("RP", mlwinlev[i], sep = ""), RP.names)
    if (length(RPx.pos) != 0) {
      cat(paste(rep("-", 50), collapse = "-"), "\n")
      RPx.names <- gsub(paste("RP+", mlwinlev[i], "+\\_", sep = ""), "", RP.names[RPx.pos])
      RPx <- as.matrix(RP.print[, RPx.pos], nrow = 4)
      printcol0 <- align2left("        ", RPx.names)
      printcol1 <- align2right("Coef.", format(round(RPx[1, ], digits), nsmall = digits))
      printcol2 <- align2right("Std. Err.", format(round(RPx[2, ], digits), nsmall = digits))
      printcol5 <- align2right("[95% Cred.", format(round(RPx[3, ], digits), nsmall = digits))
      printcol6 <- align2right("Interval]", format(round(RPx[4, ], digits), nsmall = digits))
      printcol7 <- align2right("ESS", format(round(RPx[5, ]), nsmall = 0))
      cat("The random part estimates at the", levID2[i], "level:", "\n")
      for (i in 1:(ncol(RPx) + 1)) {
        cat(printcol0[i], " ", printcol1[i], " ", printcol2[i], " ", printcol5[i], " ", printcol6[i], " ", 
            printcol7[i], "\n")
      }
    }
  }
  cat(paste(rep("-", 50), collapse = "*"), "\n")
}

#' @rdname summary-methods-mcmc
#' @export 
setMethod("print", "mlwinfitMCMC", printMCMC)

#' @rdname summary-methods-mcmc
#' @export 
setMethod("show", signature(object = "mlwinfitMCMC"), function(object) printMCMC(object))

#' Update "mlwinfitMCMC" objects
#' @param object a valid \code{mlwinfitMCMC} class object with an R function call component named \code{call}, the expression used to create itself.
#' @param Formula. changes to the formula. This is a two sided formula where "." is substituted for existing components in the \code{Formula} component of \code{object$call}.
#' @param levID. changes to the specifications of level ID(s).
#' @param estoptions. changes to the specifications of a list of options used for estimating the model.
#' @param ...  additional arguments to the call, or arguments with changed values.
#' @param keep.order a logical value indicating whether the terms should keep their positions.
#' @param evaluate  if \code{TRUE} (the default) the new call is evaluated;
#' otherwise the call is returned as an unevaluated expression.
#' @return either a new updated \code{mlwinfitMCMC} class object, or else an unevaluated expression for creating such an object.
#' @export 
setMethod("update", signature(object = "mlwinfitMCMC"), updateMLwiN)

#' "mlwinfitMCMC" model formula
#' @param x See \code{\link[stats]{formula}}
#' @param env See \code{\link[stats]{formula}}
#' @param ... Other arguments; see \code{\link[stats]{formula}}
#' @export 
setMethod("formula", "mlwinfitMCMC", function(x, env = parent.frame(), ...) {
  as.formula(x@Formula)
})

#' Extract the coefficient vector from "mlwinfitMCMC" objects
#' @param object An \code{\link{mlwinfitMCMC-class}} object
#' @param ... Other arguments
#' @seealso \code{\link[stats]{coef}}
#' @export 
setMethod("coef", signature(object = "mlwinfitMCMC"), function(object, ...) {
  c(object@FP, object@RP)
})

#' @rdname coef-mlwinfitMCMC-method
#' @export 
setMethod("coefficients", signature(object = "mlwinfitMCMC"), function(object, ...) {
  coef(object)
})

#' Extract the approximate variance-covariance matrix from "mlwinfitMCMC" objects
#' @param object An \code{\link{mlwinfitMCMC-class}} object
#' @param ... Other arguments
#' @seealso \code{\link[stats]{vcov}}
#' @export 
setMethod("vcov", signature(object = "mlwinfitMCMC"), function(object, ...) {
  m <- matrix(0, nrow(object@FP.cov) + nrow(object@RP.cov), ncol(object@FP.cov) + ncol(object@RP.cov))
  colnames(m) <- c(colnames(object@FP.cov), colnames(object@RP.cov))
  rownames(m) <- c(rownames(object@FP.cov), rownames(object@RP.cov))
  m[colnames(object@FP.cov), rownames(object@FP.cov)] <- object@FP.cov
  m[colnames(object@RP.cov), rownames(object@RP.cov)] <- object@RP.cov
  m
})

#' Returns the fitted values from "mlwinfitMCMC" objects.
#' @param object An \code{\link{mlwinfitMCMC-class}} object.
#' @param ... Other arguments
#' @seealso \code{\link[stats]{fitted}}
#' @export 
setMethod("fitted", signature(object = "mlwinfitMCMC"), function(object, ...) {
  predict(object, type = "response")
})

#' @rdname fitted-mlwinfitMCMC-method
#' @export 
setMethod("fitted.values", signature(object = "mlwinfitMCMC"), function(object, ...) {
  fitted(object)
})

#' Returns the residual data from "mlwinfitMCMC" objects.
#' @param object An \code{\link{mlwinfitMCMC-class}} object
#' @param ... Other arguments.
#' @seealso \code{\link[stats]{residuals}}
#' @export 
setMethod("residuals", signature(object = "mlwinfitMCMC"), function(object, ...) {
  form <- Formula.translate(object@Formula, object@D, object@data)
  if (!is.list(form$resp)) {
    D <- object@D
    indata <- object@data
    tval <- fitted(object)
    if (D[1] == "Binomial") {
      tval <- tval * indata[, D[3]]
    }
    if (D[1] == "Poisson") {
      if (!is.na(D[3])) {
        tval <- tval + indata[, D[3]]
      }
    }
    object@data[[form$resp]] - tval
  } else {
    warning("residuals only implemented for univariate models")
    NULL
  }
})

#' @rdname residuals-mlwinfitMCMC-method
#' @export 
setMethod("resid", signature(object = "mlwinfitMCMC"), function(object, ...) {
  residuals(object)
})

#' Returns the predicted data from "mlwinfitMCMC" objects.
#' @param object An \code{\link{mlwinfitMCMC-class}} object.
#' @param newdata data frame for which to evaluate predictions
#' @param params a character vector specifying the parameters to use in evaluating predictions. If \code{NULL}, \code{names(object[["FP"]])} is used by default.
#' @param type when this has the value \code{"link"} (default) the linear predictor is returned. 
#' When \code{type="terms"} each component of the linear predictor is returned seperately. When \code{type="response"} predictions on the scale of the response are returned.
#' @param se.fit logical. When this is \code{TRUE} (not default) standard error estimates are returned for each prediction.
#' @param terms if \code{type="terms"}, which terms (default is all terms), a character vector.
#' @param ... Other arguments
#' @seealso \code{\link[stats]{predict}}
#' @export 
setMethod("predict", signature(object = "mlwinfitMCMC"), function(object, newdata = NULL, params = NULL, type = "link", se.fit = FALSE, 
                                              terms = NULL, ...) {
  if (is.null(newdata)) {
    indata <- object@data
  } else {
    indata <- Formula.translate(object@Formula, object@D, object@data)$indata
  }
  if (is.null(params)) {
    fp.names <- names(FP <- object@FP)
  } else {
    fp.names <- params
  }
  x.names <- sub("FP_", "", fp.names)
  if (type == "link") {
    tval <- as.vector(as.matrix(indata[x.names]) %*% as.matrix(object@FP[fp.names])[, 1])
    if (se.fit) {
      # seval <- as.vector(sqrt(diag(as.matrix(indata[x.names]) %*% as.matrix(object@FP.cov[fp.names, fp.names]) %*%
      # t(as.matrix(indata[x.names])))))
      seval <- as.vector(sqrt(rowSums(as.matrix(indata[x.names]) %*% as.matrix(object@FP.cov[fp.names, fp.names]) * 
                                        indata[x.names])))
      return(list(fit = tval, se.fit = seval))
    } else {
      return(tval)
    }
  }
  if (type == "terms") {
    if (!is.null(terms)) {
      x.names <- terms
      fp.names <- paste0("FP_", terms)
    }
    tval <- as.matrix(t(t(indata[x.names]) * object@FP[fp.names]))
    if (se.fit) {
      seval <- as.matrix(sqrt(t(t(indata[x.names]^2) * diag(object@FP.cov[fp.names, fp.names]))))
      return(list(fit = tval, se.fit = seval))
    } else {
      return(tval)
    }
  }
  if (type == "response") {
    tval <- as.vector(as.matrix(indata[x.names]) %*% as.matrix(object@FP[fp.names]))
    D <- object@D
    if (D[1] == "Normal" || D[1] == "Multivariate Normal") {
      return(tval)
    }
    if (D[1] == "Binomial") {
      if (D[2] == "logit") {
        antilogit <- function(x) {
          exp(x)/(1 + exp(x))
        }
        return(antilogit(tval) * indata[, D[3]])
      }
      if (D[2] == "probit") {
        return(pnorm(tval) * indata[, D[3]])
      }
      if (D[2] == "cloglog") {
        anticloglog <- function(x) {
          1 - exp(-exp(x))
        }
        return(anticloglog(tval) * indata[, D[3]])
      }
    }
    if (D[1] == "Poisson") {
      if (is.na(D[3])) {
        return(exp(tval))
      } else {
        return(exp(tval + indata[, D[3]]))
      }
    }
    if (D[1] == "Negbinom") {
      if (is.na(D[3])) {
        return(exp(tval))
      } else {
        return(exp(tval + indata[, D[3]]))
      }
    }
    if (D[1] == "Mixed") {
    }
    if (D[1] == "Multinomial") {
    }
    warning("link function transformation not yet implemented")
    return(NULL)
  }
})

#' Returns the number of used observations from "mlwinfitMCMC" objects.
#' @param object An \code{\link{mlwinfitMCMC-class}} object.
#' @param ... Other arguments.
#' @seealso \code{\link[stats]{nobs}}
#' @export 
setMethod("nobs", signature(object = "mlwinfitMCMC"), function(object, ...) {
  object@Nobs
}) 
