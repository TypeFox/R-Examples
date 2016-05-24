.onAttach <- function(lib, pkg) {
  #checkGCDversion()
  packageStartupMessage("This is paleofire v",utils::packageDescription("paleofire",field="Version"),appendLF = TRUE)
} 





#' Check GCD package install
#' 
#' Check if GCD package is installed and up to date to ensure always using the
#' most up to date GCD version. devtools package is required: on Windows
#' install Rtools.exe depending on your R version
#' http://cran.r-project.org/bin/windows/Rtools/
#' 
#' Last GCD database version is donwloaded and installed using:
#' 
#' \code{library(devtools)}
#' 
#' \code{install_github("GCD",username="paleofire",ref="master")}
#' 
#' @author O. Blarquez
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' \dontrun{checkGCDversion()}
#' 
checkGCDversion <- function() {
  # Check to see if installed
  if (!"GCD" %in% utils::installed.packages()[, 1]) {
    Checks <- "Failed"
  } else {
    oldones=old.packages()
    if ("GCD" %in% oldones[,1]) {
      Checks <- "Failed"
    } else Checks <- "Passed"
  }
  switch(
    Checks,
    Passed = { message("Everything looks OK! GCD up to date: v", packageVersion("GCD")) },
    Failed = {
#       ans = readline(
#         "GCD is either outdated or not installed. Update now? (y/n) ")
#       if (ans != "y")
#         return(invisible())
      packageStartupMessage("GCD is either outdated or not installed. Installing...")
      install.packages("GCD")
    })
  # packageStartupMessage("This is paleofire v",utils::packageDescription("paleofire",field="Version"),appendLF = TRUE)
}


loessGCV <- function (x) {
  ## Modified from code by Michael Friendly
  ## http://tolstoy.newcastle.edu.au/R/help/05/11/15899.html
  if (!(inherits(x,"loess"))) stop("Error: argument must be a loess object")
  ## extract values from loess object
  span <- x$pars$span
  n <- x$n
  traceL <- x$trace.hat
  sigma2 <- sum(resid(x)^2) / (n-1)
  gcv  <- n*sigma2 / (n-traceL)^2
  result <- list(span=span, gcv=gcv)
  result
}

bestLoess <- function(model, spans = c(.05, .95)) {
  f <- function(span) {
    mod <- update(model, span = span)
    loessGCV(mod)[["gcv"]]
  }
  result <- optimize(f, spans)
  result
}

