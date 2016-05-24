#' User Friendly Evaluation of Dose-Effect Experiments using Litchfield-Wilcoxon
#' and Probit Methods
#'
#' User friendly evaluation of dose-effect experiments using automated
#'   Litchfield Wilcoxon (1949) and probit estimation methods.
#'   This function has been tailored for non-R users with input data set up in a
#'   particular way (see Details).
#' @param rawfile
#'   A character scalar specifying the path of the input data as a csv file.
#'     If NULL, default, the user will be prompted to browse to a file using a
#'     menu.
#' @param descrcolz
#'   A numeric vector, the column numbers to use as the description of the
#'     test, default 1:4.
#' @param saveplots
#'   A logical scalar indicating if plots should be saved to a pdf file,
#'     default TRUE.  See details.
#' @param showplots
#'   A logical scalar indicating if plots should be shown on screen,
#'     default FALSE.  See details.
#' @param saveresults
#'   A logical scalar indicating if results should be saved to a csv file,
#'     default TRUE.
#'     The csv file is given the same name (plus the suffix "Smry") and is
#'     placed in the same directory as the input file.
#' @param showresults
#'   A logical scalar indicating if results should be printed to the console,
#'   default TRUE.
#'   These results include the chi-squared statistic, degrees of freedom, and
#'   p-value for the Litchfield Wilcoxon method.
#' @param returnresults
#'   A logical scalar indicating if results should be returned by the function,
#'     default FALSE.
#' @return
#'   If \code{returnresults=TRUE}, a data frame with 11 rows per test and
#'     2 more columns than the input data. Three columns from the input data are
#'     not included (\code{TFM Conc. (mg/L)}, \code{No. Tested}, and
#'     \code{No. Dead}). Five columns are added: the parameter (\code{param}),
#'     the method used (\code{method}), the estimate (\code{estimate}), and the
#'     95\% confidence interval of the estimate (\code{lower95ci} and
#'     \code{upper95ci})
#' @import
#'   graphics grDevices tcltk utils
#' @export
#' @references
#'   Litchfield, JT Jr. and F Wilcoxon.  1949.
#'     A simplified method of evaluating dose-effect experiments.
#'     Journal of Pharmacology and Experimental Therapeutics 96(2):99-113.
#'     \href{http://jpet.aspetjournals.org/content/96/2/99.abstract}{[link]}.
#' @details
#'   The input data must include at least these seven columns,
#'     with these names in the header row:
#'   \itemize{
#'     \item \code{Test ID} = A character or numeric vector,
#'       the unique identifier for each test
#'     \item \code{Source} = A character vector, the source of the chemical
#'     \item \code{Batch} = A character or numeric vector, the batch of the
#'       chemical
#'     \item \code{Species} = A character vector, the species tested
#'     \item \code{TFM Conc. (mg/L)} = A numeric vector, the concentration of
#'       TFM in mg/L
#'     \item \code{No. Tested} = A numeric vector, the number of animals tested
#'     \item \code{No. Dead} = A numeric vector, the number of animals dead
#'   }
#'
#'   The input data are expected to be summarized by dose.
#'     If duplicate doses are provided, an error will be thrown.
#'
#'   Both \code{saveplots} and \code{showplots} may be set to FALSE, in which
#'     case no plots will be produced.  But both may not be set to TRUE.
#' @examples
#' \dontrun{
#' LWP()
#' }

LWP <- function(rawfile=NULL, descrcolz=1:4, saveplots=TRUE, showplots=FALSE,
  saveresults=TRUE, showresults=TRUE, returnresults=FALSE) {

  if (sum(abs(as.integer(descrcolz) - descrcolz))>0) {
    stop("descrcolz must be a vector of integers.")
  }
  if (!is.logical(saveplots) | length(saveplots)!=1) {
    stop("saveplots must be a logical scalar.")
  }
  if (!is.logical(showplots) | length(showplots)!=1) {
    stop("showplots must be a logical scalar.")
  }
  if (!is.logical(saveresults) | length(saveresults)!=1) {
    stop("saveresults must be a logical scalar.")
  }
  if (!is.logical(showresults) | length(showresults)!=1) {
    stop("showresults must be a logical scalar.")
  }
  if (!is.logical(returnresults) | length(returnresults)!=1) {
    stop("returnresults must be a logical scalar.")
  }
  if (saveplots & showplots) {
    stop("Either saveplots or showplots (or both) must be FALSE.")
  }

  oldopt <- as.logical(options("stringsAsFactors"))
  options(stringsAsFactors = FALSE)

  ### bring in the data ###
  # allow user to choose raw data file from menu
  if (is.null(rawfile)) rawfile <- tk_choose.files(default="*.csv", multi=FALSE)
  # read in the data, fill in the blanks
  rawdat <- read.csv(rawfile, as.is=TRUE)
  rawdat2 <- data.frame(lapply(rawdat, fill))
  rawcolz <- match(c("TFM.Conc...mg.L.", "No..Tested", "No..Dead"),
    names(rawdat2))

  if (any(descrcolz<1) | any(descrcolz>dim(rawdat)[2])) {
    stop("descrcolz should contain integers between 1 and ",
      dim(rawdat)[2], ".")
  }
  if (any(is.na(match(c("Test.ID", "Source", "Batch", "Species", "Date", "pH",
    "Temp..C.", "DO..mg.L.", "TFM.Conc...mg.L.", "No..Tested", "No..Dead"),
    names(rawdat))))) {
    stop("Input data must include all of these variables: Test ID, Source,
      Batch, Species, TFM Conc. (mg/L), No. Tested, No. Dead.")
  }
  dose.nona <- rawdat2$TFM.Conc...mg.L.[!is.na(rawdat2$TFM.Conc...mg.L.)]
  if (sum(duplicated(dose.nona))>0) {
    stop("TFM Conc. (mg/L) should be a vector of unique values, with no duplicates")
  }

  # use the input filename to name the output files
  filesegs <- strsplit(rawfile, "/")[[1]]
  L <- length(filesegs)
  filename <- filesegs[L]
  if(L>1) {
    dirname <- paste0(paste(filesegs[-L], collapse="/"), "/")
  } else {
    dirname <- ""
  }
  prefix <- strsplit(filename, ".csv")[[1]]
  smryname <- paste0(prefix, "Smry.csv")

  if (saveplots) {
    pdfname <- paste0(prefix, "Smry.pdf")
    pdf(file = paste0(dirname, pdfname), paper="letter")
    }

  ### fit LW and probit models to the data

  # unique test IDs
  sut <- sort(unique(rawdat2$Test.ID))

  # empty list in which to put results
  results <- vector("list", length(sut))

  for(i in seq(along=sut)) {
    df <- rawdat2[rawdat2$Test.ID==sut[i], ]
    descr <- paste(df[1, descrcolz], collapse=", ")
    mydat <- with(df, dataprep(dose=TFM.Conc...mg.L., ntot=No..Tested,
      nfx=No..Dead))

    fLW <- LWestimate(fitLWauto(mydat), mydat)
    fp <- fitprobit(mydat)
    pctalive <- c(25, 50, 99.9)

    pm <- predlinear(pctalive, fLW)
    estimate <- c(fLW$params, pm[, "ED"], fLW$LWest["S"])
    param <- c(names(estimate[1:2]), paste0("ED", pctalive), "S")
    method <- rep("Auto Litchfield-Wilcoxon", length(param))
    lower95ci <- c(NA, NA, pm[1, "lower"], fLW$LWest["lower"], pm[3, "lower"],
      fLW$LWest["lowerS"])
    upper95ci <- c(NA, NA, pm[1, "upper"], fLW$LWest["upper"], pm[3, "upper"],
      fLW$LWest["upperS"])
    smryLW <- data.frame(param, method, estimate, lower95ci, upper95ci)

    Pr <- do.call(rbind, lapply(pctalive, predprobit, fp))
    cp <- coefprobit(fp)
    row.names(Pr) <- paste0("ED", pctalive)
    fpc <- fp$coef
    if (!fp$converged) fpc[1:2] <- NA
    estimate <- c(fpc, Pr[, "ED"])
    param <- names(estimate)
    method <- rep("Probit", length(param))
    lower95ci <- c(cp["ilower"], cp["slower"], Pr[, "lower"])
    upper95ci <- c(cp["iupper"], cp["supper"], Pr[, "upper"])
    smryPr <- data.frame(param, method, estimate, lower95ci, upper95ci)

    smry <- rbind(smryLW, smryPr)
    n <- dim(smry)[1]
    results[[i]] <- cbind(df[rep(1, n), -rawcolz], smry)

    if (showresults) {
      # print the results to the screen
      cat("\n\n\n")
      cat(paste0("Test ", i, ":   ", descr, "\n"))
      cat("\nLitchfield Wicoxon method\n\n")
      print(fLW$chi$chi)
      cat("\n")
      print(format(smryLW[, -2], 2, nsmall=2, digits=0), row.names=FALSE)
      cat("\nProbit method\n\n")
      print(format(smryPr[, -2], 2, nsmall=2, digits=0), row.names=FALSE)
      }

    if (saveplots | showplots) {
      # plot the results to a pdf file
      par(mar=c(4, 4, 2, 1))
      plotDELP(mydat, main=descr, ylab="Mortality  (%)")
      if (!is.na(fpc[1])) abline(fpc, lty=2, col="red")
      abline(fLW$params)
      # notes on graph
      right <- 0.8 * (par("usr")[2] - par("usr")[1]) + par("usr")[1]
      lwsel <- substring(smry$method, 1, 1)=="A"
      rows <- match(c("ED25", "ED50", "ED99.9", "S"), smry$param[lwsel])
      text(right, -1.2, "Litchfield Wilcoxon", font=2)
      text(right, -seq(1.6, 2.8, 0.4), c("ED25", "ED50", "ED99.9", "LW Slope"),
        adj=1)
      text(right, -seq(1.6, 2.8, 0.4), paste("  ",
        formatC(smry$estimate[lwsel][rows], digits=3, flag="#")), adj=0)
      left <- 0.2 * (par("usr")[2] - par("usr")[1]) + par("usr")[1]
      psel <- substring(smry$method, 1, 1)=="P"
      rows <- match(c("ED25", "ED50", "ED99.9"), smry$param[psel])
      text(left, 2.9, "Probit  (dashed)", font=2, col="red")
      text(left, seq(2.5, 1.7, -0.4), c("ED25", "ED50", "ED99.9"), adj=1,
        col="red")
      text(left, seq(2.5, 1.7, -0.4),
        paste("  ", formatC(smry$estimate[psel][rows], digits=3, flag="#")),
        adj=0, col="red")
      }
    }
  if (saveplots) graphics.off()
  # save the results to a csv file
  smrydat <- do.call(rbind, results)
  if (saveresults) write.csv(smrydat, paste0(dirname, smryname),
    row.names=FALSE)
  if (returnresults) {
    return(smrydat)
  } else {
    invisible()
  }
  if (showresults) {
    # print a header to the screen
    cat("\n\n\n")
    cat("Rounded results are printed to the screen for convenience.\n")
    cat("No need to copy or print them though, because they are saved in:\n")
    cat("     ", paste0(dirname, smryname), "\n")
    cat('Note that "S" is the slope defined by',
      "Litchfield and Wilcoxon (1949).\n\n")
    }
  options(stringsAsFactors = oldopt)
  }
