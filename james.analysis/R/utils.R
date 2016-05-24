#' Open PDF file
#' 
#' Open a PDF file in the standard PDF viewer.
#' 
#' Opens the given PDF \code{file} in the default PDF viewer. Under Windows, 
#' this is done by calling \code{shell.exec} and whatever program associated 
#' with the PDF file extension will be used. On Unix (including Mac OS X) the 
#' function will use the program named in the option "pdfviewer" (see 
#' \code{help(options)} for information on how this option is set).
#' 
#' @param file character vector: relative path to the PDF file to be opened
#' @noRd
openPDF <- function(file) 
{
  # check input
  if (is.null(file) || is.na(file) || !is.character(file)){
    stop("'file' should be of type character vector")
  }
  # quote file name for in case it contains spaces
  file <- paste('"', file, '"', sep="")
  # implementation depends on operating system
  OST <- .Platform$OS.type
  if (OST == "windows") {
    shell.exec(file)
  } else if (OST == "unix") {
    # check whether pdfviewer is set
    pdf <- getOption("pdfviewer")
    msg <- NULL
    if (is.null(pdf)) {
      msg <- 'getOption("pdfviewer") is NULL'
    } else if (is.na(pdf)) {
      msg <- 'getOption("pdfviewer") is NA'
    } else if (!is.character(pdf)) {
      msg <- 'getOption("pdfviewer") should be of type "character"'
    } else if (nchar(pdf) == 0) {
      msg <- 'getOption("pdfviewer") is ""'
    }
    if (!is.null(msg)) {
      stop(msg, "; please set pdf viewer using 'options(pdfviewer = ...)'")
    }
    # create command to open pdf
    cmd <- paste(pdf, file)
    # run command
    system(cmd)
  } else {
    stop("Unknown operating system (not Windows nor Unix).")
  }
}

# Wrapper to output formatted strings.
printf <- function(...){
  cat(sprintf(...))
}

# Wrapper to check if a vector is sorted
is.sorted <- function(v){
  return(!is.unsorted(v))
}

####### DATA HANDLING #######

# Get the name of the single problem contained in the given data.
# If the data contains results for more than one problem, an error
# is thrown.
getSingleProblem <- function(data){
  all.problems <- getProblems(data)
  if(length(all.problems) > 1){
    stop("'data' contains more than one problem; please specifiy 'problem'")
  } else {
    # return name of only problem
    return(all.problems[1])
  }
}

# Get the name of the single search applied to the given problem
# in the given data. If the data contains results for more than
# one search for this problem, an error is thrown. Argument
# 'problem' can be omitted if 'data' contains results for a
# single problem only.
getSingleSearch <- function(data, problem){
  all.searches <- getSearches(data, problem)
  if(length(all.searches) > 1){
    msg <- sprintf("'data' contains more than one search for problem \"%s\"; please specify 'search'", problem)
    stop(msg)
  } else {
    # return name of only search
    return(all.searches[1])
  }
}

# Check if values are being maximized or minimized for a given problem,
# by comparing the first and last obtained value. If contradicting results
# are found for different searches or runs solving the same problem, an
# error is thrown. Argument 'problem' can be omitted if 'data' contains
# results for a single problem only.
isMinimizing <- function(data, problem){
  # fall back to single problem if missing
  if(missing(problem)){
    problem <- getSingleProblem(data)
  }
  # get searches
  searches <- getSearches(data, problem)
  # check maximization/minimization for each run of each search
  minimizing <- c()
  for(search in searches){
    runs <- getSearchRuns(data, problem, search)
    for(run in runs){
      first.value <- head(run$values, 1)
      last.value <- tail(run$values, 1)
      if(first.value > last.value){
        minimizing <- c(minimizing, TRUE)
      } else if (first.value < last.value){
        minimizing <- c(minimizing, FALSE)
      }
    }
  }
  # process results
  if(length(minimizing) == 0){
    # all convergence curves are "flat"; can be both (and does not matter)
    # default to maximizing
    return(FALSE)
  } else {
    minimizing <- unique(minimizing)
    if(length(minimizing) != 1){
      msg <- sprintf(paste("contradicting results for problem \"%s\";",
                           "some search runs were maximizing,",
                           "other minimizing"),
                     problem)
      stop(msg)
    } else {
      return(minimizing[1])
    }
  }
}

### TIME UNITS ###

# get abbreviation of a time unit
getTimeUnitAbbreviation <- function(time.unit = c("milliseconds", "seconds", "minutes", "hours")){
  abbr <- switch(match.arg(time.unit),
                 "milliseconds" = "ms",
                 "seconds" = "sec",
                 "minutes" = "min",
                 "hours" = "h")
  return(abbr)
}
# get the number of milliseconds contained in a specific time unit
getTimeUnitMillis <- function(time.unit = c("milliseconds", "seconds", "minutes", "hours")){
  ms <- switch(match.arg(time.unit),
                 "milliseconds" = 1,
                 "seconds" = 1000,
                 "minutes" = 60000,
                 "hours" = 3600000)
  return(ms)
}

### RELEASE QUESTIONS ###

release_questions <- function() {
  c(
    "Has the version number been updated? Check DESCRIPTION, package docs and NEWS.md.",
    "Has the release date been set? Check DESCRIPTION, package docs and NEWS.md",
    "Do the package docs refer to the correct version(s) of the extensions module?",
    "Are the changes listed in NEWS.md up to date?",
    "Pushed everything to GitHub?"
  )
}





