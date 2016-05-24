asFunction <- function(expr, envir=parent.frame(), ...) {
  eval(substitute(function() x, list(x=expr)), envir=envir, ...)
}

#' @importFrom utils installed.packages
findBasePkgs <- local({
  pkgs <- NULL
  function() {
    if (length(pkgs) > 0L) return(pkgs)
    data <- installed.packages()
    isBase <- (data[,"Priority"] %in% "base")
    pkgs <<- rownames(data)[isBase]
    pkgs
  }
})

isBasePkgs <- function(pkgs) {
  pkgs <- gsub("^package:", "", pkgs)
  pkgs %in% findBasePkgs()
}

# cf. is.primitive()
is.base <- function(x) {
  if (typeof(x) != "closure") return(FALSE)
  isBasePkgs(environmentName(environment(x)))
}

# cf. is.primitive()
is.internal <- function(x) {
  if (typeof(x) != "closure") return(FALSE)
  body <- deparse(body(x))
  any(grepl(".Internal", body, fixed=TRUE))
}

## Emulates R internal findVar1mode() function
## https://svn.r-project.org/R/trunk/src/main/envir.c
where <- function(x, where=-1, envir=if (missing(frame)) { if (where < 0) parent.frame(-where) else as.environment(where) } else sys.frame(frame), frame, mode="any", inherits=TRUE) {
  tt <- 1
  ## Validate arguments
  stopifnot(is.environment(envir))
  stopifnot(is.character(mode), length(mode) == 1L)
  inherits <- as.logical(inherits)
  stopifnot(inherits %in% c(FALSE, TRUE))

  ## Search
  while (!identical(envir, emptyenv())) {
    if (exists(x, envir=envir, mode=mode, inherits=FALSE)) return(envir)
    if (!inherits) return(NULL)
    envir <- parent.env(envir)
  }

  NULL
}


## From R.utils 2.0.2 (2015-05-23)
hpaste <- function(..., sep="", collapse=", ", lastCollapse=NULL, maxHead=if (missing(lastCollapse)) 3 else Inf, maxTail=if (is.finite(maxHead)) 1 else Inf, abbreviate="...") {
  if (is.null(lastCollapse)) lastCollapse <- collapse

  # Build vector 'x'
  x <- paste(..., sep=sep)
  n <- length(x)

  # Nothing todo?
  if (n == 0) return(x)
  if (is.null(collapse)) return(x)

  # Abbreviate?
  if (n > maxHead + maxTail + 1) {
    head <- x[seq(length=maxHead)]
    tail <- rev(rev(x)[seq(length=maxTail)])
    x <- c(head, abbreviate, tail)
    n <- length(x)
  }

  if (!is.null(collapse) && n > 1) {
    if (lastCollapse == collapse) {
      x <- paste(x, collapse=collapse)
    } else {
      xT <- paste(x[1:(n-1)], collapse=collapse)
      x <- paste(xT, x[n], sep=lastCollapse)
    }
  }

  x
} # hpaste()


## From future 0.11.0
trim <- function(s) {
  sub("[\t\n\f\r ]+$", "", sub("^[\t\n\f\r ]+", "", s))
} # trim()


## From future 0.11.0
hexpr <- function(expr, trim=TRUE, collapse="; ", maxHead=6L, maxTail=3L, ...) {
  code <- deparse(expr)
  if (trim) code <- trim(code)
  hpaste(code, collapse=collapse, maxHead=maxHead, maxTail=maxTail, ...)
} # hexpr()
