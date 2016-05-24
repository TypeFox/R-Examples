.distrTEstoptions <- list(
                      MaxNumberofPlottedEvaluationDims = 6,
                      MaxNumberofPlottedEvaluations = 6,
                      MaxNumberofSummarizedEvaluations = 15,
                      MaxNumberofPrintedEvaluations = 15)
  
distrTEstoptions <- function(...) {
  if (nargs() == 0) return(.distrTEstoptions)
  current <- .distrTEstoptions
  temp <- list(...)
  if (length(temp) == 1 && is.null(names(temp))) {
    arg <- temp[[1]]
    switch(mode(arg),
           list = temp <- arg,
           character = return(.distrTEstoptions[arg]),
           stop("invalid argument: ", sQuote(arg)))
  }
  if (length(temp) == 0) return(current)
  n <- names(temp)
  if (is.null(n)) stop("options must be given by name")
  changed <- current[n]
  current[n] <- temp
  if (sys.parent() == 0) env <- asNamespace("distrTEst") 
  else                   env <- parent.frame()
  assign(".distrTEstoptions", current, envir = env)
  invisible(current)
}

getdistrTEstOption <- function(x) distrTEstoptions(x)[[1]]
