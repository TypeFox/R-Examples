"steppath" <-
function (y, ..., max.blocks) UseMethod("steppath")

"steppath.default" <-
function (y, x = 1:length(y), x0 = 2 * x[1] - x[2], max.cand = NULL, family = c("gauss", "gaussvar", "poisson", "binomial", "gaussKern"), param = NULL, weights = rep(1, length(y)), cand.radius = 0, ..., max.blocks = max.cand)
{
  family <- match.arg(family)
  cand <- stepcand(y = y, x = x, x0 = x0, max.cand = max.cand, family = family, param = param, weights = weights, cand.radius = cand.radius)
  steppath(cand, max.blocks = min(sum(!is.na(cand$number)), max.blocks))
}

"steppath.stepcand" <-
function (y, ..., max.blocks = sum(!is.na(y$number)))
{
  cand <- y
  # selected candidates
  sel.cand <- which(!is.na(y$number))
  if(max.blocks > length(sel.cand)) stop("number of blocks max.blocks may not exceed number of candidates")
  if(max.blocks < 1) stop("number of blocks max.blocks must be positive")
  algo <- switch(attr(y, "family"),
#     gaussInhibitBoth = "gaussInhibit",
    attr(y, "family")
  )
  ret <- switch(algo,
    gauss = {
      .Call(.pathGauss, y$cumSum, y$cumSumSq, y$cumSumWe, max.blocks)
    },
    gaussvar = {
      .Call(.pathGaussVar, y$cumSumSq, y$cumSumWe, max.blocks)
    },
    gaussKern = {
#       .Call(.pathGaussInhibit, y$cumSum, y$cumSumSq, y$cumSumWe, max.blocks, as.integer(attr(y, "param")$inhibit["start"]), as.integer(attr(y, "param")$inhibit["middle"]), as.integer(attr(y, "param")$inhibit["end"]))
      y <- y[sel.cand,]
      crows <- 1:nrow(cand)
      before <- attr(y, "param")$inhibit["middle"]
      after <- attr(y, "param")$inhibit["middle"]
      p <- .Call(.pathGaussCut, y$bcumSum, y$bcumSumSq, y$bcumSumWe, y$acumSum, y$acumSumSq, y$acumSumWe, max.blocks, as.integer(before), as.integer(after))
      p$path <- lapply(p$path, function(ri) {
        ri <- sel.cand[ri] # turn into rows of y
        k <- length(ri)
        # further improve jump selection
        if(k > 1) for(i in 1:(k - 1)) {
          if(i == 1) {
            sel <- which(crows <= ri[i + 1]) # between neighbouring jumps
            ri[i] <- sel[.Call(.forwardGaussInhibit, cand$cumSum[sel], cand$cumSumSq[sel], cand$cumSumWe[sel], as.integer(2), as.integer(before), as.integer(before), as.integer(after))$rightIndex[1]]
          } else {
            sel <- which(crows > ri[i - 1] & crows <= ri[i + 1]) # between neighbouring jumps
            sel.ri <- .Call(.forwardGaussInhibit, cand$cumSum[sel] - cand$cumSum[ri[i - 1]], cand$cumSumSq[sel] - cand$cumSumSq[ri[i - 1]], cand$cumSumWe[sel] - cand$cumSumWe[ri[i - 1]], as.integer(2), as.integer(before), as.integer(before), as.integer(after))$rightIndex
            if(length(sel.ri) == 2) ri[i] <- sel[sel.ri[1]]
          }
        }
        ri
      })
      p
    },
#     gaussInhibit = {
#       .Call(.pathGaussInhibit, y$cumSum, y$cumSumSq, y$cumSumWe, max.blocks, as.integer(attr(y, "param")["start"]), as.integer(attr(y, "param")["middle"]), as.integer(attr(y, "param")["end"]))
#     },
    poisson = {
      .Call(.pathPoisson, y$cumSum, y$cumSumWe, max.blocks)
    },
    binomial = {
      .Call(.pathBinom, as.integer(attr(y, "param")), y$cumSum, y$cumSumWe, max.blocks)
    },
    stop("unknown family")
  )
  ret$cand <- cand
  class(ret) <- c("steppath", class(ret))
  ret
}

"[[.steppath" <-
function(x, i)
{
  if(is.character(i)) {
    x[[i]]
  } else {
    ri <- x$path[[i]]
    ret <- stepfit(cost = x$cost[i], family = attr(x$cand, "family"), value = rep(NA, length(ri)), param = attr(x$cand, "param"), 
      leftEnd = x$cand$leftEnd[c(0, ri[-length(ri)]) + 1], rightEnd = x$cand$rightEnd[ri], x0 = attr(x$cand, "x0"),
      leftIndex = x$cand$leftIndex[c(0, ri[-length(ri)]) + 1], rightIndex = x$cand$rightIndex[ri])
    algo <- switch(attr(x$cand, "family"),
#       gaussInhibit = "gauss",
#       gaussInhibitBoth = "gauss",
      attr(x$cand, "family")
    )
    switch(algo, 
      gauss = {
      ret$cumSum <- x$cand$cumSum[ri]
      ret$cumSumSq <- x$cand$cumSumSq[ri]
      ret$cumSumWe <- x$cand$cumSumWe[ri]
      ret$value <- diff(c(0, ret$cumSum)) / diff(c(attr(ret, "x0"), ret$rightIndex))
      },
      gaussKern = {
        param <- attr(x$cand, "param")
        if(i > 1) {
          r <- x$cand$rightIndex[ri]
          ir <- r[-length(r)] # inner positions
          iri <- ri[-length(ri)]
          n <- r[length(r)] # last position
          kl <- length(param$kern)
          kj <- param$jump
          s <- param$step
          d <- c(ir - kj, n) - c(0, ir + kl - kj) + c(rep(sum((1 - s)^2), length(r) - 1), 0) + c(0, rep(sum(s^2), length(r) - 1))
          ld <- length(d)
          sd <- sum( s * ( 1 - s ) )
          Xy <- c(0, x$cand$lXy[iri]) + x$cand$rcXy[ri] - c(0, x$cand$lcXy[iri]) + x$cand$rXy[ri]
          # X'X : diagonal carries inner length of block + sum(s^2) + sum((1-s)^2), secondary diagonal is constant sum( s ( 1 - s ) )
#         if(attr(x$cand, "Matrix")) { # use sparsity
#           require(Matrix, quietly = TRUE)
#           XX <- bandSparse(ld, ld, 0:1, list(d, rep(sd, ld - 1)), symmetric = T)
#         } else {
          XX <- diag(d)
          XX[cbind(1:(ld-1), 2:ld)] <- sd
          XX[cbind(2:ld, 1:(ld-1))] <- sd
#         }
          ret$value <- as.numeric(solve(XX, Xy))
        } else {
          ret$value <- x$cand$cumSum[ri] / ri
        }
      },
      poisson = {
        ret$cumSum <- x$cand$cumSum[ri]
        ret$cumSumWe <- x$cand$cumSumWe[ri]
        ret$value <- diff(c(0, ret$cumSum)) / diff(c(attr(ret, "x0"), ret$rightIndex))
      },
      binomial = {
        ret$cumSum <- x$cand$cumSum[ri]
        ret$cumSumWe <- x$cand$cumSumWe[ri]
        ret$value <- diff(c(0, ret$cumSum)) / diff(c(attr(ret, "x0"), ret$rightIndex)) / attr(ret, "param")
      },
      stop("unknown family")
    )
    ret
  }
}

"length.steppath" <-
function(x)
{
  length(x$path)
}

"print.steppath" <-
function(x, ...)
{
  cat("\n")
  cat("Fitted solution path of step functions\n\n")
  cat("max. number of blocks:", length(x), "\n")
  cat("family:", attr(x$cand, "family"), "\n")
  if(!is.null(attr(x$cand, "param"))) {
    cat("parameter:\n")
    print(attr(x$cand, "param"))
  }
  cat("\n")
}

"logLik.steppath" <-
function(object, df = NULL, nobs = object$cand$rightIndex[nrow(object$cand)], ...)
{
  if(is.null(df)) df <- 1:length(object$path)
  ret <- switch(attr(object$cand, "family"),
    gauss = {
      param <- attr(object$cand, "param")
      if(is.null(param)) {
	df <- df + 1 # variance has also been estimated
	-nobs / 2 * (1 + log(2 * pi * object$cost / nobs))
      } else {
	-nobs / 2 * log(2 * pi * param^2) - object$cost / 2 / param^2
      }
    },
    gaussKern = stop("family gaussKern not implemented"),
    -object$cost
  )
  attr(ret, "df") <- df
  attr(ret, "nobs") <- nobs
  class(ret) <- "logLik"
  ret
}
