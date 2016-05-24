#-------------------------------------------------------------------------------
# File: boa.core.q
# Description: Library of core functions for the Bayesian Output Analysis
#    Program (BOA)
# Author: Brian J. Smith <brian-j-smith@uiowa.edu>
#-------------------------------------------------------------------------------

.boa.chain <- NULL
.boa.par <- NULL
.boa.pardesc <- NULL
.boa.version <- NULL

boa.chain <- function(...)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   if(nargs() == 0)  return(.boa.chain)
   temp <- list(...)
   if((length(temp) == 1) && is.null(names(temp))) {
      arg <- temp[[1]]
      switch(mode(arg),
         list = temp <- arg,
         character = return(.boa.chain[[arg]]),
         temp <- NULL
      )
   }
   idx <- intersect(names(temp), names(.boa.chain))
   if(length(idx) > 0) {
      current <- .boa.chain
      current[idx] <- temp[idx]
      assignInNamespace(".boa.chain", current, "boa")
   } else {
      cat("Warning: invalid arguments\n")
   }
   invisible()
}


boa.init <- function(recover = FALSE)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   recover <- recover && !is.null(.boa.par) && !is.null(.boa.pardesc) &&
	           !is.null(.boa.chain) && !is.null(.boa.version)

   if(!recover) {   
      assignInNamespace(".boa.par", list(
		   acf.lags     = c(1, 5, 10, 50),
         alpha        = 0.05,
         bandwidth    = function(x)  0.5 * diff(range(x)) / (log(length(x)) + 1),
         batch.size   = 50,
         dev.list     = numeric(0),
         ASCIIext     = ".txt",
         gandr.bins   = 20,
         gandr.win    = 0.5,
         geweke.bins  = 10,
         geweke.first = 0.1,
         geweke.last  = 0.5,
         handw.error  = 0.1,
         kernel       = "gaussian",
         legend       = TRUE,
         path         = "",
         par          = list(),
         plot.mfdim   = c(3, 2),
         plot.new     = FALSE,
         plot.onelink = FALSE,
         quantiles    = c(0.025, 0.5, 0.975),
         randl.error  = 0.005,
         randl.delta  = 0.001,
         randl.q      = 0.025,
         title        = TRUE), "boa")
      assignInNamespace(".boa.pardesc", structure(rbind(
         c("Convergence", "Brooks, Gelman, & Rubin", "Alpha Level", "alpha", ""),
         c("Convergence", "Brooks, Gelman, & Rubin", "Window Fraction", "gandr.win", ""),
         c("Convergence", "Geweke", "Window 1 Fraction", "geweke.first", ""),
         c("Convergence", "Geweke", "Window 2 Fraction", "geweke.last", ""),
         c("Convergence", "Heidelberger & Welch", "Accuracy", "handw.error", ""),
         c("Convergence", "Heidelberger & Welch", "Alpha Level", "alpha", ""),
         c("Convergence", "Raftery & Lewis", "Accuracy", "randl.error", ""),
         c("Convergence", "Raftery & Lewis", "Alpha Level", "alpha", ""),
         c("Convergence", "Raftery & Lewis", "Delta", "randl.delta", ""),
         c("Convergence", "Raftery & Lewis", "Quantile", "randl.q", ""),
         c("Convergence Plot", "Brooks & Gelman", "Number of Bins", "gandr.bins", ""),
         c("Convergence Plot", "Brooks & Gelman", "Window Fraction", "gandr.win", ""),
         c("Convergence Plot", "Gelman & Rubin", "Alpha Level", "alpha", ""),
         c("Convergence Plot", "Gelman & Rubin", "Number of Bins", "gandr.bins", ""),
         c("Convergence Plot", "Gelman & Rubin", "Window Fraction", "gandr.win", ""),
         c("Convergence Plot", "Geweke", "Alpha Level", "alpha", ""),
         c("Convergence Plot", "Geweke", "Number of Bins", "geweke.bins", ""),
         c("Convergence Plot", "Geweke", "Window 1 Fraction", "geweke.first", ""),
         c("Convergence Plot", "Geweke", "Window 2 Fraction", "geweke.last", ""),
         c("Descriptive", "Statistics", "ACF Lags", "acf.lags", ""),
         c("Descriptive", "Statistics", "Alpha Level", "alpha", ""),
         c("Descriptive", "Statistics", "Batch Size", "batch.size", ""),
         c("Descriptive", "Statistics", "Quantiles", "quantiles", ""),
         c("Descriptive Plot", "Density", "Bandwidth", "bandwidth",
           "This defines the standard deviation of the smoothing kernel"),
         c("Descriptive Plot", "Density", "Kernel", "kernel",
           "Possible kernels are gaussian, rectangular, triangular, or cosine"),
         c("Import", "Files", "Working Directory", "path",
           "Use forward slashes (\'/\') as directory separators and omit a terminating slash"),
         c("Import", "Files", "ASCII File Ext", "ASCIIext", ""),
         c("Plot", "Graphics", "Legend", "legend", "Include plot legends (T/F)"),
         c("Plot", "Graphics", "Title", "title", "Include plot title (T/F)"),
         c("Plot", "Graphics", "Keep Previous Plots", "plot.new", ""),
         c("Plot", "Graphics", "Plot Layout", "plot.mfdim",
           paste("A vector of the form \'c(nr, nc)\' giving the number of\n",
                 "rows and columns to include in each plot", sep="")),
         c("Plot", "Graphics", "Plot Chains Separately", "plot.onelink", ""),
         c("Plot", "Graphics", "Graphical Parameters", "par", "")),
         dimnames = list(NULL, c("group", "method", "desc", "par", "note"))),
			"boa")
      assignInNamespace(".boa.chain", list(
		   master         = list(),
         master.support = list(),
         work           = list(),
         work.support   = list(),
         work.sync      = TRUE), "boa")
      assignInNamespace(".boa.version", list(
		   name     = "BOA",
         major    = 1,
         minor    = 1,
         revision = 8,
         system   = version$system), "boa")
      boa.license()
   }
   
   invisible()
}


boa.load <- function(name, envir = globalenv())
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   loaded <- FALSE
   if(exists(name, envir = envir)) {
      obj <- get(name, envir = envir)
      if(is.character(obj$version$name) && is.numeric(obj$version$major) &&
         is.numeric(obj$version$minor)) {
         ver <- boa.version()
         loaded <- (obj$version$name == ver$name) &&
                   ((obj$version$major < ver$major) ||
                       ((obj$version$major == ver$major)
                       && (obj$version$minor >= ver$minor)))
      }
      if(loaded) {
         boa.par(obj$par)
         boa.par(dev.list = numeric(0))
         boa.chain(obj$chain)
      } else {
         cat("Warning: object is incompatible with this version of BOA\n")
      }
   } else {
      cat("Warning: object not found\n")
   }

   return(loaded)
}


boa.par <- function(...)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   if(nargs() == 0)  return(.boa.par)
   temp <- list(...)
   if((length(temp) == 1) && is.null(names(temp))) {
      arg <- temp[[1]]
      switch(mode(arg),
         list = temp <- arg,
         character = return(.boa.par[[arg]]),
         temp <- NULL
      )
   }
   changed <- NULL
   globals <- names(temp)
   idx <- is.element(globals, names(.boa.par))
   if(!all(idx)) {
      cat("Warning: invalid arguments\n")
      print(globals[!idx])
   }
   if(any(idx)) {
      globals <- globals[idx]
      pclass <- unlist(lapply(.boa.par[globals], "data.class"))
      idx <- unlist(lapply(temp[globals], "data.class")) == pclass
      if(!all(idx)) {
         cat("Warning: arguments must be of type\n")
         print(pclass[!idx])
      }
      if(any(idx)) {
         globals <- globals[idx]
         current <- .boa.par
         changed <- current[globals]
         current[globals] <- temp[globals]
         assignInNamespace(".boa.par", current, "boa")
      }
   }
   invisible(changed)
}


boa.quit <- function()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   cat("Exiting BOA session...\n")
   assignInNamespace(".boa.par", NULL, "boa")
   assignInNamespace(".boa.pardesc", NULL, "boa")
   assignInNamespace(".boa.chain", NULL, "boa")
   assignInNamespace(".boa.version", NULL, "boa")
   invisible()
}


boa.save <- function(name, envir = globalenv(), replace = FALSE)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
	if (replace || !is.element(name, objects(envir))) {
      assign(name, list(par = boa.par(), chain = boa.chain(),
                        version = boa.version()), envir = envir)
		saved <- TRUE
	} else {
		saved <- FALSE
	}
	
	saved
}
