#' Multiple comparison procedures to the means of a factor using
#' the studentized range and midrange distributions.
#'
#' \code{MRtest} applies the Skott-Knott midrange, Skott-Knott range,
#'     Student-Newman-Keuls midrange and Tukey midrange tests. These are new
#'     tests for multiple comparisons proposed by the authors (2015), that are in publication
#'     fase.
#'
#' @param y Model (aov or lm), numeric vector containing the
#'     response variable or the mean of the treatments.
#' @param trt Constant (y = model) or a vector containing
#'     the treatments.
#' @param dferror Degrees of freedom of the Mean Square Error.
#' @param mserror Mean Square Error.
#' @param replication Number de repetitions of the treatments in the experiment.
#' For unbalanced data should be informed the harmonic mean of repetitions.
#' This argument should be informed only if \code{ismean = TRUE}.
#' @param alpha Significant level. The default is \code{alpha = 0.05}.
#' @param main Title of the analysis.
#' @param MCP Allows choosing the multiple comparison test;
#'     the \emph{defaut} is "all". This option will go perform all tests.
#'     However, the options are: the Skott-Knott midrange test ("SKM"),
#'     the Skott-Knott Range test ("SKR"), the Student-Newman-Keuls
#'     midrange test ("SNKM") and the Tukey midrange test ("TM").
#' @param ismean Logic. If \code{FALSE} (default), the \code{y} argument represents
#' a model (aov or lm) or a numeric vector containing the response variable. If \code{TRUE} the \code{y} argument represents the
#' mean of treatments.
#' @return \code{MRtest} returns the print of a list of results. First,
#'     the summary of \code{y}. Second, the statistics
#'     of the test chosen. And finally, the  mean group results for each test.
#'     If \code{MRtest} function is stored
#'     in an object, the results will be printed and
#'     also stored in the object.
#' @details The \code{MCP} argument allows you to choose various tests
#'     of multiple comparisons at once. For example,
#'     \code{MCP = c("SKM", "SKR")}, and so on.
#' @examples
#' # Simulated data (completely randomized design)
#'
#' # Response variable
#' rv <- c(100.08, 105.66, 97.64, 100.11, 102.60, 121.29, 100.80,
#'         99.11, 104.43, 122.18, 119.49, 124.37, 123.19, 134.16,
#'         125.67, 128.88, 148.07, 134.27, 151.53, 127.31)
#'
#' # Treatments
#' treat <- factor(rep(LETTERS[1:5], each = 4))
#'
#' # Anova
#' res     <- anova(aov(rv~treat))
#' DFerror <- res$Df[2]
#' MSerror <- res$`Mean Sq`[2]
#'
#' # Loading the midrangeMCP package
#' library(midrangeMCP)
#'
#' # applying the tests
#' results <- MRtest(y = rv,
#'                   trt = treat,
#'                   dferror = DFerror,
#'                   mserror = MSerror,
#'                   alpha = 0.05,
#'                   main = "Multiple Comparison Procedure: SKM test",
#'                   MCP = c("SKM"))
#'
#' # Other option for the MCP argument is "all". All tests are used.
#'
#' results$Groups     # Results of the tests
#' results$Statistics # Main arguments of the tests
#' results$Summary    # Summary of the response variable
#'
#' # Using the y argument as aov or lm model
#' res  <- aov(rv~treat)
#'
#' MRtest(y = res, trt = "treat", alpha = 0.05,
#'        main = "Multiple Comparison Procedure: SKM test",
#'        MCP = c("SKM"))
#'
#' # For unbalanced data: It will be used the harmonic mean of
#' #                       the number of experiment replicates
#'
#' # Using the previous example
#' rv <- rv[-1]
#' treat <- treat[-1]
#'
#' res  <- lm(rv~treat) # Linear model
#'
#' # Multiple comparison procedure: SKR test
#' MRtest(y = res, trt = "treat", alpha = 0.05,
#'        main = "Multiple Comparison Procedure: SKR test",
#'        MCP = c("SKR"))
#'
#' # Assuming that the available data are the averages
#' #  of the treatments and the analysis of variance
#'
#' # Analysis of Variance Table
#'
#' # Response: rv
#' #            Df Sum Sq Mean Sq F value    Pr(>F)
#' # treat      4 4135.2 1033.80  14.669 4.562e-05 ***
#' # Residuals 15 1057.1   70.47
#'
#' mean.treat <- c(100.87, 105.95, 117.62, 127.97, 140.30)
#' treat <- factor(LETTERS[1:5])
#' DFerror <- 15
#' MSerror <- 70.47488
#' replic <- 4
#'
#' MRtest(y = mean.treat,
#'        trt = treat,
#'        dferror = DFerror,
#'        mserror = MSerror,
#'        replication = replic,
#'        alpha = 0.05,
#'        main = "Multiple Comparison Procedure: SKM test",
#'        MCP = c("SKM"),
#'        ismean = TRUE)
#'
#' @import "stats" "utils" "graphics" "SMR"
#' @export
MRtest <- function(y, trt = NULL, dferror = NULL, mserror = NULL, replication = NULL, alpha = 0.05, main = NULL,
                   MCP = "all", ismean = FALSE){
  #####################################################
  #Defensive programming
  if (ismean == TRUE) {
    if (is.null(replication)) {
      stop("The replication argument must be informed", call. = FALSE)
    }
    if (!is.numeric(y)) {
      stop("The y argument must be numeric", call. = FALSE)
    }
    if (!is.factor(trt)) {
      stop("The trt argument must be factor", call. = FALSE)
    }
  }
  if (is.numeric(y)) {
    if (!is.numeric(y)) {
      stop("The y argument must be numeric", call. = FALSE)
    }
    if (!is.factor(trt)) {
      stop("The trt argument must be factor", call. = FALSE)
    }
    if (length(y) != length(trt)) {
      stop("The y and trt arguments must have same length", call. = FALSE)
    }
    if (is.null(dferror)) {
      stop("The dferror argument must be informed", call. = FALSE)
    }
    if (is.null(mserror)) {
      stop("The mserror argument must be informed", call. = FALSE)
    }
  }
  #####################################################
  if (all(MCP == "all")){
    MCP = c("SKM", "SKR", "SNKM", "TM")
  }
  #####################################################
  #Defensive programming
  if(is.null(trt)){
    stop("The trt argument is required", call. = FALSE)
  }
  mcps <- c("SKM", "SKR", "SNKM", "TM")
  nas  <- pmatch(MCP, mcps)
  if (any(is.na(nas))) {
    stop("The options for the MCP argument are 'SKM',
         'SKR', 'SNKM' and 'TM'", call. = FALSE)
  }
  ################################################
  name.y   <- paste(deparse(substitute(y)))
  name.trt <- paste(deparse(substitute(trt)))
  if (is.null(main)) {
    main  <- paste(name.y, "~", name.trt)
  }
  #y: "aov" or "lm" class
  if ("aov" %in% class(y) | "lm" %in% class(y)) {
    if (is.null(main)) {
      main <- y$call
    }
    A <- y$model
    dferror <- df.residual(y)
    mserror <- deviance(y)/dferror
    y <- A[, 1]
    if (length(trt) > length(names(A)[-1]) + 1) {
      stop("The length of the trt argument is invalid", call. = FALSE)
    }
    ipch <- pmatch(trt, names(A))
    if (any(is.na(ipch))) {
      ncolnam <- length(colnames(A))
      optrt   <- colnames(A)[2]
      if (ncolnam > 2) {
        for (i in 3:ncolnam) {
          optrt <- paste(optrt, names(A)[i], sep = " ")
        }
      }
      stop("The options for the trt argument are\n Options: ", optrt, "\n Note: Observe which of the options have signified practice", call. = FALSE)
    }
    ipch  <- pmatch(trt, names(A))
    nipch <- length(ipch)
    ntrt  <- length(names(A))
    if (any(is.na(ipch))) {
      optrt <- names(A)[2]
      for(i in 3:ntrt) {
        optrt <- paste(optrt, names(A)[i], sep = " ")
      }
      stop("Any of the options of trt argument is wrong \n Options: ", optrt, call. = FALSE)
    }
    name.t <- names(A)[ipch][1]
    trt    <- A[, ipch]
    if (nipch > 1) {
      trt <- A[, ipch[1]]
      for (i in 2:nipch) {
        name.t <- paste(name.t, names(A)[ipch][i], sep = ":")
        trt <- paste(trt, A[, ipch[i]], sep = ":")
      }
    }
    name.y <- names(A)[1]
  }
  cat("MCP's based on distributions of the studentized midrange and range\n\n")
  cat("Study: ", main, "\n\n")

  if (ismean == TRUE) {
    data.na    <- subset(data.frame(y, trt), is.na(y) == FALSE)
    means      <- round(y, 2)
    n          <- length(y)
    rn         <- rep(replication, n)
    rh         <- 1/mean(1/rn)
    Mean       <- mean(y)
    CV         <- sqrt(mserror) * 100/Mean
  } else {
    data.na    <- subset(data.frame(y, trt), is.na(y) == FALSE)
    means      <- round(tapply(data.na[, 1], data.na[, 2],  "mean"), 2)
    std        <- round(tapply(data.na[, 1], data.na[, 2], "sd"), 2)
    rn         <- tapply(data.na[, 1], data.na[, 2], "length")
    mi         <- round(tapply(data.na[, 1], data.na[, 2], "min"), 2)
    ma         <- round(tapply(data.na[, 1], data.na[, 2], "max"), 2)
    n          <- nrow(means)
    rh         <- 1/mean(1/rn)
    Mean       <- mean(data.na[, 1])
    CV         <- sqrt(mserror) * 100/Mean
  }


  ################################################
  #Observation for unbalanced data
  if (length(unique(rn)) != 1) {
    cat("Unbalanced data: It will be used the harmonic mean of
        \n the number of experiment replicates \n")
  }
  ################################################

  if (ismean == TRUE) {
    summarydata <- data.frame(Means  = means,
                              r      = rn)
    cat("Summary:\n")
    print(summarydata)
  } else {
    summarydata <- data.frame(Means  = means,
                              std    = std,
                              r      = rn,
                              Min    = mi,
                              Max    = ma)
    cat("Summary:\n")
    print(summarydata)
  }

  if (length(unique(rn)) != 1) {
    cat("\nHarmonic mean of the number of experiment replicates ",
        rh, "\n")
  }
  #DMS midrange
  if(any(MCP == "SNKM")) {
    nnn      <- 2:n
    aaa      <- rep(1 - alpha / 2, times = (n - 1))
    MRq.snkm <- SMR::qSMR(aaa, nnn, dferror) # Studentized midrange of the SNKM test
    MRq      <- MRq.snkm[n - 1] # Studentized midrange of the SKM/TM test
    dms.int  <- MRq * sqrt(mserror / rh) # Internal DMS of the SKM/TM test
    dms      <- dms.int + sqrt(0.5) * sqrt(mserror / rh) / n^0.5 # add of mean square error/nmed^0.5
    dms      <- c(dms, dms.int) # External DMS of the SKM/TM test
    dmsint   <- MRq.snkm * sqrt(mserror / rh) # Internal DMS of the SNKM test
    dms.snkm <- dmsint + sqrt(0.5) * sqrt(mserror / rh) / n^0.5 # External DMS of the SNKM test
  } else{
    MRq     <- SMR::qSMR(1 - alpha / 2, n, dferror) # Studentized midrange of the SKM/TM test
    dms.int <- MRq * sqrt(mserror / rh) # Internal DMS of the SKM/TM test
    dms     <- dms.int + sqrt(0.5) * sqrt(mserror/rh) / n^0.5 # External DMS of the SKM/TM test
    dms     <- c(dms, dms.int) # DMS of the SKM/TM test
  }

  #DMS range:
  Rq   <- qtukey(1 - alpha, n, dferror) # Studentized range of the SKR test
  dms.range <- Rq * sqrt(mserror / rh) # DMS of the SKR test

  #Initial statistics:
  statistics.SKM  <- NA
  statistics.SKR  <- NA
  statistics.SNKM <- NA
  statistics.TM   <- NA

  #Initial groups:
  test.SKM  <- NA
  test.SKR  <- NA
  test.SNKM <- NA
  test.TM   <- NA


  # Skott-Knott Midrange Test
  if (any(MCP == "SKM")){
    cat("\nSkott-Knott Midrange Test\n\n")
    statistics <- data.frame(Exp.Mean = Mean,
                             CV      = CV,
                             MSerror = mserror,
                             Df      = dferror,
                             n       = n,
                             Stud.Midrange = MRq,
                             Ext.DMS = dms[1],
                             Int.DMS = dms[2])
    cat("Statistics: \n")
    rownames(statistics) <- " "
    statistics.SKM <- statistics
    print(statistics)
    test <- SKMtest(y, trt, n, dferror, mserror, alpha, dms)
    test[1] <- round(test[1], 2)
    test.SKM <- test
    cat("\nGroups: \n")
    print(test)
  }

  # Skott-Knott Range Test
  if (any(MCP == "SKR")){
    cat("\nSkott-Knott Range Test\n\n")
    statistics <- data.frame(Exp.Mean = Mean,
                             CV      = CV,
                             MSerror = mserror,
                             Df      = dferror,
                             n       = n,
                             Stud.Range = Rq,
                             DMS = dms.range)
    cat("Statistics: \n")
    rownames(statistics) <- " "
    statistics.SKR   <- statistics
    print(statistics)
    test <- SKRtest(y, trt, n, dferror, mserror, alpha, dms.range)
    test[1] <- round(test[1], 2)
    test.SKR <- test
    cat("\nGroups: \n")
    print(test)
  }

  # Student-Newman-Keuls (SNK) Midrange Test
  if (any(MCP == "SNKM")){
    cat("\nSNK Midrange Test\n\n")
    statistics <- data.frame(Exp.Mean = Mean,
                             CV      = CV,
                             MSerror = mserror,
                             Df      = dferror,
                             n       = nnn,
                             Stud.Midrange = MRq.snkm,
                             DMS = dms.snkm)
    statistics <- statistics[order(statistics[, 6], decreasing = FALSE),]
    cont <- seq(1, n - 1, by = 1)
    rownames(statistics) <- cont
    rownames(statistics) <- rownames(MRq.snkm, do.NULL = FALSE, prefix = "comp")
    cat("Statistics: \n")
    statistics.SNKM <- round(statistics, 4)
    print(round(statistics, 4))
    test <- SNKMtest(y, trt, n, dferror, mserror, alpha, dms.snkm)
    test[1] <- round(test[1], 2)
    test.SNKM <- test
    cat("\nGroups: \n")
    print(test)
  }

  # Tukey Midrange Test
  if (any(MCP == "TM")){
    cat("\nTukey Midrange Test\n\n")
    statistics <- data.frame(Exp.Mean = Mean,
                             CV      = CV,
                             MSerror = mserror,
                             Df      = dferror,
                             n       = n,
                             Stud.Midrange = MRq,
                             Ext.DMS = dms[1],
                             Int.DMS = dms[2])
    cat("Statistics: \n")
    rownames(statistics) <- " "
    statistics.TM <- statistics
    print(statistics)

    test <- TMtest(y, trt, n, dferror, mserror, alpha, dms)
    test[1] <- round(test[1], 2)
    test.TM <- test
    cat("\nGroups: \n")
    print(test)
  }

  #All statistics
  stat.tests <- list(Statistics.SKM  = statistics.SKM,
                     Statistics.SKR  = statistics.SKR,
                     Statistics.SNKM = statistics.SNKM,
                     Statistics.TM   = statistics.TM)

  #All groups
  group.tests <- list(group.SKM  = test.SKM,
                      group.SKR  = test.SKR,
                      group.SNKM = test.SNKM,
                      group.TM   = test.TM)
  ################
  # Output results
  ################

  alltest    <- 1:4
  nas        <- nas[order(nas, na.last = NA)]
  ntest      <- alltest[- nas]

  if (length(nas) == 4) {
    statistics <- stat.tests
  } else {
    statistics <- stat.tests[- ntest]
  }

  if(length(nas) == 4) {
    grouptest  <- group.tests
  } else{
    grouptest  <- group.tests[- ntest]
  }

  output <- list(Summary = summarydata,
                    Groups  = grouptest,
                 Statistics = statistics,
                     Tests  = MCP)
  invisible(output)
}
