"stepfit" <-
function(cost, family, value, param = NULL, leftEnd, rightEnd, x0, leftIndex = leftEnd, rightIndex = rightEnd)
{
  ret <- stepblock(value, leftEnd, rightEnd, x0)
  ret$leftIndex <- leftIndex
  ret$rightIndex <- rightIndex
  attr(ret, "cost") <- cost
  attr(ret, "family") <- family
  attr(ret, "param") <- param
  class(ret) <- c("stepfit", class(ret))
  ret
}

"print.stepfit" <-
function(x, ...)
{
  cat("\n")
  cat("Fitted step function of family", attr(x, "family"), "containing", nrow(x), "blocks\n\n")
  cat("domain: (", attr(x, "x0"), ",", x$rightEnd[nrow(x)], "]\n")
  cat("range:  [", min(x$value), ",", max(x$value), "]\n")
  cat("cost:", attr(x, "cost"), "\n")
  if(!is.null(attr(x, "param"))) {
    cat("parameter:\n")
    print(attr(x, "param"))
  }
  cat("\n")
}

"[.stepfit" <- 
function (x, i, j, drop = if(missing(i)) TRUE else if(missing(j)) FALSE else length(j) == 1, refit = FALSE) 
{
  ret <- NextMethod("[.")
  # refit
  if(!identical(refit, FALSE)) {
    if(missing(i)) i <- 1:nrow(x)
    ret$value <- switch(attr(x, "family"),
      gauss = diff(c(0, ret$cumSum)) / diff(c(attr(ret, "x0"), ret$rightIndex)),
      gaussKern = if(nrow(ret) == 1 ) diff(c(0, ret$cumSum)) / diff(c(attr(ret, "x0"), ret$rightIndex)) else {
        param <- attr(x, "param")
        r <- ret$rightIndex
        ir <- r[-length(r)] # inner positions
        n <- r[length(r)] # last position
        kl <- length(param$kern)
        kj <- param$jump
        s <- param$step
        d <- c(ir - kj, n) - c(0, ir + kl - kj) + c(rep(sum((1 - s)^2), length(r) - 1), 0) + c(0, rep(sum(s^2), length(r) - 1))
        ld <- length(d)
        sd <- sum( s * ( 1 - s ) )
        XX <- diag(d)
        XX[cbind(1:(ld-1), 2:ld)] <- sd
        XX[cbind(2:ld, 1:(ld-1))] <- sd
        dif <- diff(c(0, ret$rightIndex))
        close <- min(dif) <= kl
        Xy <- c(0, ret$lXy[-nrow(ret)]) + ret$rcXy - c(0, ret$lcXy[-nrow(ret)]) + ret$rXy
        if(identical(refit, TRUE) | !close) {
          if(close) warning("jumps closer than filter length")
        } else {
          ss <- c(0, s, 1)
          revss <- c(1, rev(s), 0)
          for(i in which(dif < kl)) {
            if(i > 1 & i < n) {
              # compute design matrix for neighbouring blocks
              neigh <- which(r + kj > r[i-1] & c(0, ir) <= r[i] + kl - kj )
              neighn <- length(neigh)
              neighi <- max(r[neigh[1]] - kj + 1, 1):min(r[neigh[neighn] - 1] + kl - kj, n)
              tX <- outer(1:neighn, neighi, function(k,ind) {
#                 print(data.frame(k = k, neighk = neigh[k], rneighk = r[neigh[k]], ind = ind, test = ind <= r[neigh[k]] - kj, s = pmin(pmax(c(0,r)[neigh[k]] + kl - kj + 1 - ind, 0), kl + 1) + 1, revss = revss[pmin(pmax(c(0,r)[neigh[k]] + kl - kj + 1 - ind, 0), kl + 1) + 1], s1 = pmin(pmax(ind + kj - r[neigh[k]], 0), kl + 1) + 1, ss = ss[pmin(pmax(ind + kj - r[neigh[k]], 0), kl + 1) + 1]))
                revss[pmin(pmax(c(0,r)[neigh[k]] + kl - kj + 1 - ind, 0), kl + 1) + 1] - ss[pmin(pmax(ind + kj - r[neigh[k]], 0), kl + 1) + 1]
              })
#               print(tX)
              iXX <- tX %*% t(tX)
              # apply to off-diagonal elements including or crossing i
              XX[neigh[neigh <= i],neigh[neigh >= i]] <- iXX[which(neigh <= i),which(neigh >= i)]
              XX[neigh[neigh >= i],neigh[neigh <= i]] <- iXX[which(neigh >= i),which(neigh <= i)]
              Xy[i] <- tX[which(neigh == i),] %*% refit[neighi]
            }
          }
        }
        ret$value <- as.numeric(solve(XX, Xy))
      },
      poisson = diff(c(0, ret$cumSum)) / diff(c(attr(ret, "x0"), ret$rightIndex)),
      binomial = diff(c(0, ret$cumSum)) / diff(c(attr(ret, "x0"), ret$rightIndex)) / attr(ret, "param")
    )
    class(ret) <- class(x)
  }
  ret
}

"plot.stepfit" <-
function(x, dataspace = TRUE, ...)
{
  if(attr(x, "family") == "binomial" && dataspace) {
    x$value <- x$value * attr(x, "param")
  }
  NextMethod()
}

"lines.stepfit" <-
function(x, dataspace = TRUE, ...)
{
  if(attr(x, "family") == "binomial" && dataspace) {
    x$value <- x$value * attr(x, "param")
  }
  NextMethod()
}

"fitted.stepfit" <-
function(object, ...)
{
  if(attr(object, "family") == "gaussKern" && nrow(object) > 1) {
    k <- attr(object, "param")$kern
    l <- length(k)
    j <- attr(object, "param")$jump
    ret <- rep(object$value, object$rightIndex - object$leftIndex + 1 + c(l - j - 1, rep(0, length(object$rightIndex) - 2), j))
    convolve(ret, rev(k), conj = TRUE, type = "filter")
  } else if(attr(object, "family") == "binomial") {
    rep(object$value * attr(object, "param"), object$rightIndex - object$leftIndex + 1)
  } else {
    rep(object$value, object$rightIndex - object$leftIndex + 1)
  }
}

"residuals.stepfit" <-
function(object, y, ...)
{
  fit <- fitted(object)
  if(length(fit) != length(y)) stop("data and fit differ in length")
  switch(attr(object, "family"),
    binomial = {
      y - attr(object, "param") * fit
    },
    gaussvar = {
      ifelse(fit == 0, ifelse(y^2 == 0, 0, Inf), log(y^2 / fit)) # sort of residuals: log( y^2 / sigma^2 )
    },
    {
      y - fit
    }
  )
}

"logLik.stepfit" <-
function(object, df = NULL, nobs = object$rightIndex[nrow(object)], ...)
{
  family <- switch(attr(object, "family"),
    gaussKern = "gauss",
#     gaussInhibit = "gauss",
#     gaussInhibitBoth = "gauss",
    attr(object, "family")
  )
  ret <- switch(family,
    gauss = {
      if(is.null(df)) df <- nrow(object) + 1 # variance has also been estimated
      -nobs / 2 * (1 + log(2 * pi * attr(object, "cost") / nobs))
    },
    {
      if(is.null(df)) df <- nrow(object)
      attr(object, "cost")
    }
  )
  attr(ret, "df") <- df
  attr(ret, "nobs") <- nobs
  class(ret) <- "logLik"
  ret
}
