# solarBaeysAvgPolygenic
#
# @export
solarBaeysAvgPolygenic <- function(traits, covlist = "1", covlist.test,
  household = as.logical(NA),
  data, dir,
  polygenic.settings = "",  polygenic.options = "",
  bayesavg.max = as.integer(0),
  verbose = 0,
  ...) 
{
  ### step 1: args
  mc <- match.call()
  
  stopifnot(!missing(traits))
  stopifnot(length(traits) == 1)
  stopifnot(!missing(covlist.test))
  stopifnot(!missing(data))
  stopifnot(class(data) == "data.frame")
  
  is.tmpdir <- missing(dir)
  
  # force `household <- FALSE` if the data set does not have the house-hold field
  if(is.na(household) & !hasHousehold(names(data))) {
    household <- FALSE
  }
  
  # set up `polygenic.settings`
  if(is.na(household)) {
    polygenic.settings <- c(polygenic.settings, "house")
  } else if(household) {
    polygenic.settings <- c(polygenic.settings, "house")
    polygenic.options <- paste(polygenic.options, "-keephouse")
  } 
  
  # check `traits`, `covlist`
  check_var_names(traits, c(covlist, covlist.test), names(data))

  out <- list(traits = traits, covlist = covlist, covlist.test = covlist.test, 
    polygenic.settings = polygenic.settings, polygenic.options = polygenic.options, 
    solar = list(model.filename = "null0.mod", phe.filename = "dat.phe", kinship = FALSE),
    call = mc, verbose = verbose)

  ### step 2: set up SOLAR dir
  if(is.tmpdir) {
    dir <- tempfile(pattern = "solarBaeysAvgPolygenic-")
  }
  if(verbose > 1) cat(" * solarBaeysAvgPolygenic: parameter `dir` is missing.\n")
  if(verbose > 1) cat("  -- temporary directory `", dir, "` used\n")
 
  df2solar(data, dir)
  
  ### step 3: run polygenic
  out <- solar_polygenic(dir, out)
  
  ### step 4: run bayes avr.
  writeLines(covlist.test, file.path(dir, "testcovars.lst"))
  
  if(length(bayesavg.max)) {
    bmax <- length(covlist.test)
  } else {
    bmax <- min(bayesavg.max, length(covlist.test))
  }

  cmd <- c(paste0("load model ", file.path(traits, out$solar$model.filename)),
    paste0("bayesavg -overwrite -cov -max ", bmax, " -list testcovars.lst"))
  
  ret <- solar(cmd, dir)
  
  # read bayes avr. results
  line <- tail(ret, 1)
  cov2 <- as.character(NA)
  pat <- "\\*\\*\\* Model with best BIC loaded:"
  line.ok <- grepl(pat, line)
  if(line.ok) {
    line <- gsub(pat, "", line)
    line <- gsub("cov", "", line)
    line <- gsub(" ", "", line)
    ind <- as.numeric(strsplit(line, "_")[[1]])
    cov2 <- try(covlist.test[ind])
  }
  covlist.bayesavg <- c(covlist, cov2)
  
  ### output
  out$covlist.bayesavg <- covlist.bayesavg
  out$solar$bayesavg$ret <- ret
  out$solar$bayesavg$files$bayesavg_cov.out <- suppressWarnings(try(readLines(file.path(dir, traits, "bayesavg_cov.out"))))

  ### clean 
  if(is.tmpdir) {
    unlink(dir, recursive = TRUE)
    if(verbose > 1) cat("  -- solarBaeysAvgPolygenic: temporary directory `", dir, "` unlinked\n")
  }
  
  oldClass(out) <- "solarBaeysAvgPolygenic"  
  return(out)
}


#--------------------
# Class
#--------------------

#' S3 class solarBaeysAvgPolygenic
#'
#' @name solarBaeysAvgPolygenicClass
#' @rdname solarBaeysAvgPolygenicClass
#'
#' @param x 
#'    An object of class \code{solarBaeysAvgPolygenic}.
#' @param object
#'    An object of class \code{solarBaeysAvgPolygenic}.
#' @param ...
#'    Additional arguments.
#' @exportClass solarBaeysAvgPolygenic

#--------------------
# Print method
#--------------------

#' @rdname solarBaeysAvgPolygenicClass
#' @export
print.solarBaeysAvgPolygenic <- function(x, ...)
{
  cat("\nCall: ")
  print(x$call)
  
  cat("\nFinal list of covariates: ")
  cat(paste(x$covlist.bayesavg, collapse = ", "), "\n")

}

#' @rdname solarBaeysAvgPolygenicClass
#' @export
summary.solarBaeysAvgPolygenic <- function(object, ...)
{
  cat("\nFile bayesavg_cov.out:\n")
  l_ply(object$solar$bayesavg$files$bayesavg_cov.out, function(x) cat(x, "\n"))
}
