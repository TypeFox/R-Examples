
se <- function(x, ...) UseMethod("se", x)

se.default <- function(x, ...) gnm::se(x, ...)

se.rc <- function(x, type=c("se", "quasi.se"), ...) {
  if(!inherits(x, "rc")) 
      stop("x must be a rc object")

  if(length(x$assoc) == 0)
      stop("x must have an association component")

  se.assoc(x$assoc, type=type, ...)
}

se.rcL <- function(x, type=c("se", "quasi.se"), ...) {
  if(!inherits(x, "rcL")) 
      stop("x must be a rcL object")

  if(length(x$assoc) == 0)
      stop("x must have an association component")

  se.assoc(x$assoc, type=type, ...)
}

se.hmskew <- function(x, type=c("se", "quasi.se"), ...) {
  if(!inherits(x, "hmskew"))
      stop("x must be a hmskew object")

  if(length(x[["assoc"]]) > 0 && length(x$assoc.hmskew) > 0)
      return(list(assoc=se.assoc(x$assoc, type=type, ...),
                  assoc.hmskew=se.assoc(x$assoc.hmskew, type=type, ...)))
  if(length(x[["assoc"]]) > 0)
      return(se.assoc(x$assoc, type=type, ...))
  else if(length(x$assoc.hmskew) > 0)
      return(se.assoc(x$assoc.hmskew, type=type, ...))
  else
      stop("x must have an association or a skew-association component")
}

se.hmskewL <- function(x, type=c("se", "quasi.se"), ...) {
  if(!inherits(x, "hmskew"))
      stop("x must be a hmskew object")

  if(length(x[["assoc"]]) > 0 && length(x$assoc.hmskew) > 0)
      return(list(assoc=se.assoc(x$assoc, type=type, ...),
                  assoc.hmskew=se.assoc(x$assoc.hmskewL, type=type, ...)))
  if(length(x[["assoc"]]) > 0)
      return(se.assoc(x$assoc, type=type, ...))
  else if(length(x$assoc.hmskew) > 0)
      return(se.assoc(x$assoc.hmskew, type=type, ...))
  else
      stop("x must have an association or a skew-association component")
}

se.yrcskew <- function(x, type=c("se", "quasi.se"), ...) {
  if(!inherits(x, "yrcskew"))
      stop("x must be a yrcskew object")

  if(length(x[["assoc"]]) > 0 && length(x$assoc.yrcskew) > 0)
      return(list(assoc=se.assoc(x$assoc, type=type, ...),
                  assoc.yrcskew=se.assoc(x$assoc.yrcskew, type=type, ...)))
  if(length(x[["assoc"]]) > 0)
      return(se.assoc(x$assoc))
  else if(length(x$assoc.yrcskew) > 0)
      return(se.assoc(x$assoc.yrcskew))
  else
      stop("x must have an association or a skew-association component")
}

se.assoc <- function(x, type=c("se", "quasi.se"), ...) {
  type <- match.arg(type)

  if(!inherits(x, "assoc"))
      stop("x must be an assoc object")

  if(x$covtype == "none" || length(x$covmat) == 0)
      stop("No covariance matrix found: use the 'se' argument when fitting model")

  if(!(ncol(x$row) == ncol(x$col) &&
       ncol(x$phi) == ncol(x$row)))
      stop("Invalid component length")

  nd <- ncol(x$phi)
  nl <- nrow(x$phi)
  nlr <- dim(x$row)[3]
  nlc <- dim(x$col)[3]
  nr <- nrow(x$row)
  nc <- nrow(x$col)

  if(nrow(x$covmat) != ncol(x$covmat) ||
     nrow(x$covmat) != nl * nd + nlr * nd * nr + nlc * nd * nc)
      stop("Covariance matrix dimensions do not match association structure")

  if(nlr != nlc)
     stop("Different number of layers for rows and columns is currently not supported")

  std.errs <- list()

  covmat <- x$covmat

  if(type == "quasi.se") {
      get.se <- function(int) qvcalc::qvcalc(covmat[int, int, drop=FALSE])$qvframe$quasiSE
  }
  else {
      get.se <- function(int) sqrt(diag(covmat[int, int, drop=FALSE]))
  }

  std.errs$phi <- x$phi
  std.errs$phi[] <- get.se(seq.int(1, nl * nd))

  std.errs$row <- x$row
  std.errs$col <- x$col

  int <- nl * nd + rep(seq.int(0, nlr * nd - 1) * (nr + nc), each=nr) + 1:nr
  std.errs$row[] <- get.se(int)

  int <- nl * nd + nr + rep(seq.int(0, nlc * nd - 1) * (nr + nc), each=nc) + 1:nc
  std.errs$col[] <- get.se(int)

#  int <- nl * nd + rep((1:nlr - 1) * (nd * nr + nd * nc), each=nd * nr) +
#                   rep((1:nd - 1) * (nr + nc), each=nr) +
#                   seq.int(1, nr) - 1

  std.errs
}

