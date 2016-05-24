#' @title jModelTest
#' @description Run jModelTest to determine appropriate substitution model.
#' 
#' @param x a list of DNA sequences.
#' @param sub.schemes number of substitution schemes to test. Can be one 
#'   of 3, 5, 7, 11, or 203.
#' @param unequal.base.freq logical. Include models with unequal base 
#'   frequencies?
#' @param prop.inv.sites logical. Include models with a proportion of 
#'   invariable sites?
#' @param rate.var number of categories for models with rate variation 
#'   among sites.
#' @param AIC,AICc,BIC,DT logical. Calculate respective information 
#'   criterion metrics?
#' @param param.imp logical. Calculate parameter importances?
#' @param model.average logical. Do model averaging and parameter importances?
#' @param numThreads Number of threads to use.
#' @param path path where \code{jModelTest.jar} is located.
#' @param java.opts options to \code{java} command line.
#' 
#' @note jModelTest is not included with \code{strataG} and must be downloaded 
#'   separately. Additionally, it must be installed such that it can be run from 
#'   the command line in the current working directory. See the vignette 
#'   for \code{external.programs} for installation instructions.
#' 
#' @references Darriba D, Taboada GL, Doallo R, Posada D. 2012. 
#'   jModelTest 2: more models, new heuristics and parallel computing. 
#'   Nature Methods 9(8), 772.\cr
#'   Guindon S and Gascuel O (2003). A simple, fast and accurate method 
#'   to estimate large phylogenies by maximum-likelihood". Systematic 
#'   Biology 52: 696-704.\cr
#'   Available at: \url{https://code.google.com/p/jmodeltest2/}
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
jmodeltest <- function(x, sub.schemes = 3, unequal.base.freq = FALSE,
  prop.inv.sites = FALSE, rate.var = NULL, AIC = FALSE, AICc = FALSE,
  BIC = FALSE, DT = FALSE, param.imp = FALSE, model.average = FALSE, 
  numThreads = 1, path = ifelse(
    .Platform$OS.type == "windows", 
    "C:/Program Files/jModelTest", 
    "/usr/local/bin/jmodeltest"
  ), java.opts = NULL) {
  
  if(inherits(x, "multidna")) {
    if(getNumLoci(x) > 1) warning("'x' is a multidna object with more than one loci. Using first locus.")
    x <- getSequences(x, loci = 1, simplify = TRUE)
  }
  
  if(!sub.schemes %in% c(3, 5, 7, 11, 203)) {
    stop("'sub.schemes' not equal to 3, 5, 7, 11, or 203.")
  }
  
  wd <- getwd()
  setwd(path)
  in.file <- write.fasta(x, file = "jModelTest.fasta")
  output.file <- gsub(".fasta", ".results.txt", in.file)
  
  modeltest.call <- paste(
    "java", java.opts, "-jar jModelTest.jar",
    "-d", in.file,
    "-s", sub.schemes,
    ifelse(unequal.base.freq, "-f", ""),
    ifelse(prop.inv.sites, "-i", ""),
    ifelse(!is.null(rate.var), paste("-g", rate.var), ""),
    ifelse(AIC, "-AIC", ""),
    ifelse(AICc, "-AICc", ""),
    ifelse(BIC, "-BIC", ""),
    ifelse(DT, "-DT", ""),
    ifelse(param.imp, "-p", ""),
    ifelse(model.average, "-v", ""),
    "-tr", numThreads,
    ">", output.file
  )
  err.code <- if(.Platform$OS.type == "unix") {
    system(modeltest.call, intern = F)
  } else {
    shell(modeltest.call, intern = F)
  }
  setwd(wd)
  
  if(err.code == 0) {
    file.remove(file.path(path, in.file))
    file.rename(file.path(path, output.file), file.path(wd, output.file))
    cat("jModelTest finished successfully\n")
    invisible(output.file)
  } else {
    warning(paste("jModelTest returned error code", err.code))
    invisible(NULL)
  }
}