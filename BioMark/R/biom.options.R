## biom.options.R: modeled after pls.options.R.
## New version August 4, 2014. (Thanks for the example BH!)

## The list of initial options, for the moment all pertaining to
## stability-based biomarkers selection

## 21-11-2011: added "lasso" as possibility for fmethods. We assume
## that the number of observations is smaller than the number of
## variables. Eventually we also take lasso to be an empty list. Let's see.

biom.options <- function(..., reset = FALSE) {
  if (reset) {
    .biom.Options$options <-
        list(max.seg = 100, oob.size = NULL, oob.fraction = .3,
             variable.fraction = .5, ntop = 10, min.present = .1,
             fmethods = c("studentt", "shrinkt", "pcr", "pls", "vip", "lasso"),
             univ.methods = c("studentt", "shrinkt"),
             lasso = list(alpha = 1, nlambda = 100),
             nset = 10000, HCalpha = .1)

    .biom.Options$options
  }

  temp <- list(...)
  if (length(temp) == 0)
      .biom.Options$options
  
  current <- .biom.Options$options
  if (length(temp) == 1 && is.null(names(temp))) {
    arg <- temp[[1]]
    switch(mode(arg),
           list = temp <- arg,
           character = return(.biom.Options$options[arg]),
           stop("invalid argument: ", sQuote(arg)))
  }
  if (length(temp) == 0) return(current)
  n <- names(temp)
  if (is.null(n)) stop("options must be given by name")
  changed <- current[n]
  current[n] <- temp
  .biom.Options$options <- current
  invisible(current)
}

.biom.Options <- new.env(parent = emptyenv())
biom.options(reset = TRUE)

