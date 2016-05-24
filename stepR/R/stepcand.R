"stepcand" <-
function(y, x = 1:length(y), x0 = 2 * x[1] - x[2], max.cand = NULL, family = c("gauss", "gaussvar", "poisson", "binomial", "gaussKern"), param = NULL, weights = rep(1, length(y)), cand.radius = 0)
{
  family <- match.arg(family)
  if(!is.null(max.cand)) if(max.cand > length(y)) stop("number of blocks max.cand may not exceed length of y")
  if(!is.null(max.cand)) if(max.cand < 1) stop("number of blocks max.cand must be positive")
  if(length(x) != length(y)) stop("x and y must of the same length")
  if(length(x0) != 1) stop("x0 must be a single numeric")
  if(length(weights) != length(y)) stop("weights and y must of the same length")
  if(family == "binomial") {
    if(is.null(param)) stop("param must be size of binomial distribution")
    if(any(y < 0) | any(y > param)) stop("only values in {0, ..., param} allowed for family binomial")
    y <- as.integer(y)
  }
  if(family == "poisson") {
    if(any(y < 0)) stop("only values in the natural numbers allowed for family poisson")
    y <- as.integer(y)
  }
  if(family == "gaussKern") {
    kl <- length(param$kern)
    # step response
    s <- param$step
    kj <- param$jump
    # inhibit blocks shorter than kernel
    param$inhibit <- c(start = kj + 1, middle = kl, end = kl - kj)
#     before <- kj + 1
#     after <- kl - kj
    before <- after <- kl
  }
#   if(family %in% c("gaussInhibit", "gaussInhibitBoth")) {
#     if(length(param) == 3 & is.null(names(param))) names(param) <- c("start", "middle", "end")
#     if(length(param) == 1) param <- c(start = param, middle = param, end = param)
#     if(any(sort(names(param)) != c("end", "middle", "start"))) stop("param for family gaussInhibit(Both) needs start, middle, and end")
#   }
  algo <- switch(family,
    gaussKern = "gaussCut",
#     gaussInhibit = "gauss",
#     gaussInhibitBoth = "gaussInhibit",
    family
  )
  switch(algo,
    gauss = {
      cumSum <- as.double(cumsum(y))
      cumSumSq <- as.double(cumsum(y^2))
      cumSumWe <- as.double(cumsum(weights))
      if(!is.null(max.cand)) {
        forward <- .Call(.forwardGauss, cumSum, cumSumSq, cumSumWe, as.integer(max.cand))
        rightIndex <- neighbours(forward$rightIndex, 1:length(y), cand.radius)
      } else {
        forward <- structure(NA, cost = NA)
        rightIndex <- 1:length(y)
      }
      ret <- stepfit(cost = attr(forward, "cost"), family = family, param = param, value = rep(NA, length(rightIndex)),
        leftEnd = c(x[1], x[rightIndex[-length(rightIndex)] + 1]), rightEnd = x[rightIndex], x0 = x0,
        leftIndex = c(1, rightIndex[-length(rightIndex)] + 1), rightIndex = rightIndex)
      ret$cumSum <- cumSum[rightIndex]
      ret$cumSumSq <- cumSumSq[rightIndex]
      ret$cumSumWe <- cumSumWe[rightIndex]
      ret$value <- diff(c(0, ret$cumSum)) / diff(c(0, ret$rightIndex))
    },
    gaussvar = {
      cumSumSq <- as.double(cumsum(y^2))
      cumSumWe <- as.double(cumsum(weights))
      if(!is.null(max.cand)) {
        forward <- .Call(.forwardGaussVar, cumSumSq, cumSumWe, as.integer(max.cand))
        rightIndex <- neighbours(forward$rightIndex, 1:length(y), cand.radius)
      } else {
        forward <- structure(NA, cost = NA)
        rightIndex <- 1:length(y)
      }
      ret <- stepfit(cost = attr(forward, "cost"), family = family, param = param, value = rep(NA, length(rightIndex)),
        leftEnd = c(x[1], x[rightIndex[-length(rightIndex)] + 1]), rightEnd = x[rightIndex], x0 = x0,
        leftIndex = c(1, rightIndex[-length(rightIndex)] + 1), rightIndex = rightIndex)
      ret$cumSumSq <- cumSumSq[rightIndex]
      ret$cumSumWe <- cumSumWe[rightIndex]
      ret$value <- diff(c(0, ret$cumSumSq)) / diff(c(0, ret$rightIndex))
    },
#     gaussInhibit = {
#       if(is.null(max.cand)) max.cand <- length(y)
#       cumSum <- as.double(cumsum(y))
#       cumSumSq <- as.double(cumsum(y^2))
#       cumSumWe <- as.double(cumsum(weights))
#       forward <- .Call(.forwardGaussInhibit, cumSum, cumSumSq, cumSumWe, as.integer(max.cand), as.integer(param["start"]), as.integer(param["middle"]), as.integer(param["end"]))
#       rightIndex <- neighbours(forward$rightIndex, 1:length(y), cand.radius)
#       ret <- stepfit(cost = attr(forward, "cost"), family = family, param = param, value = rep(NA, length(rightIndex)),
#         leftEnd = c(x[1], x[rightIndex[-length(rightIndex)] + 1]), rightEnd = x[rightIndex], x0 = x0,
#         leftIndex = c(1, rightIndex[-length(rightIndex)] + 1), rightIndex = rightIndex)
#       ret$cumSum <- cumSum[rightIndex]
#       ret$cumSumSq <- cumSumSq[rightIndex]
#       ret$cumSumWe <- cumSumWe[rightIndex]
#       ret$value <- diff(c(0, ret$cumSum)) / diff(c(0, ret$rightIndex))
#     },
    gaussCut = {
      if(is.null(max.cand)) max.cand <- length(y)
      n <- length(y)
      cumSum <- as.double(cumsum(y))
      cumSumSq <- as.double(cumsum(y^2))
      cumSumWe <- as.double(cumsum(weights))
      bcumSum <- as.double(c(rep(NA, before), cumSum[-(n + 1 - 1:before)]))
      bcumSumSq <- as.double(c(rep(NA, before), cumSumSq[-(n + 1 - 1:before)]))
      bcumSumWe <- as.double(c(rep(NA, before), cumSumWe[-(n + 1 - 1:before)]))
      acumSum <- as.double(c(cumSum[-(1:after)], rep(NA, after)))
      acumSumSq <- as.double(c(cumSumSq[-(1:after)], rep(NA, after)))
      acumSumWe <- as.double(c(cumSumWe[-(1:after)], rep(NA, after)))
      forward <- .Call(.forwardGaussCut, bcumSum, bcumSumSq, bcumSumWe, acumSum, acumSumSq, acumSumWe, as.integer(max.cand), as.integer(before), as.integer(after))
      rightIndex <- neighbours(forward$rightIndex, 1:length(y), cand.radius)
#       if(family == "gaussKern") {
#         rightIndex <- c(rightIndex[rightIndex > kj & rightIndex <= n - 1 - kl + kj], n)
#       }
      ret <- stepfit(cost = attr(forward, "cost"), family = family, param = param, value = rep(NA, length(rightIndex)),
        leftEnd = c(x[1], x[rightIndex[-length(rightIndex)] + 1]), rightEnd = x[rightIndex], x0 = x0,
        leftIndex = c(1, rightIndex[-length(rightIndex)] + 1), rightIndex = rightIndex)
      ret$cumSum <- cumSum[rightIndex]
      ret$cumSumSq <- cumSumSq[rightIndex]
      ret$cumSumWe <- cumSumWe[rightIndex]
      ret$acumSum <- acumSum[rightIndex]
      ret$acumSumSq <- acumSumSq[rightIndex]
      ret$acumSumWe <- acumSumWe[rightIndex]
      ret$bcumSum <- bcumSum[rightIndex]
      ret$bcumSumSq <- bcumSumSq[rightIndex]
      ret$bcumSumWe <- bcumSumWe[rightIndex]
      ret$value <- diff(c(0, ret$cumSum)) / diff(c(0, ret$rightIndex))
    },
    poisson =  {
      cumSum <- as.integer(cumsum(y))
      cumSumWe <- as.double(cumsum(weights))
      if(!is.null(max.cand)) {
        forward <- .Call(.forwardPoisson, cumSum, cumSumWe, as.integer(max.cand))
        rightIndex <- neighbours(forward$rightIndex, 1:length(y), cand.radius)
      } else {
        forward <- structure(NA, cost = NA)
        rightIndex <- 1:length(y)
      }
      ret <- stepfit(cost = attr(forward, "cost"), family = family, param = param, value = rep(NA, length(rightIndex)),
        leftEnd = c(x[1], x[rightIndex[-length(rightIndex)] + 1]), rightEnd = x[rightIndex], x0 = x0,
        leftIndex = c(1, rightIndex[-length(rightIndex)] + 1), rightIndex = rightIndex)
      ret$cumSum <- cumSum[rightIndex]
      ret$cumSumWe <- cumSumWe[rightIndex]
      ret$value <- diff(c(0, ret$cumSum)) / diff(c(0, ret$rightIndex))
    },
    binomial =  {
      cumSum <- as.integer(cumsum(y))
      cumSumWe <- as.double(cumsum(weights))
      if(!is.null(max.cand)) {
        forward <- .Call(.forwardBinom, as.integer(param), cumSum, cumSumWe, as.integer(max.cand))
        rightIndex <- neighbours(forward$rightIndex, 1:length(y), cand.radius)
      } else {
        forward <- structure(NA, cost = NA)
        rightIndex <- 1:length(y)
      }
      ret <- stepfit(cost = attr(forward, "cost"), family = family, param = param, value = rep(NA, length(rightIndex)),
        leftEnd = c(x[1], x[rightIndex[-length(rightIndex)] + 1]), rightEnd = x[rightIndex], x0 = x0,
        leftIndex = c(1, rightIndex[-length(rightIndex)] + 1), rightIndex = rightIndex)
      ret$cumSum <- cumSum[rightIndex]
      ret$cumSumWe <- cumSumWe[rightIndex]
      ret$value <- diff(c(0, ret$cumSum)) / diff(c(0, ret$rightIndex)) / param
    },
    stop("unknown family")
  )
  if(!is.null(max.cand)) {
    ret$number <- rep(NA, nrow(ret))
    ret$number[rightIndex %in% forward$rightIndex] <- forward$number
  } else {
    ret$number <- rep(0, nrow(ret))
  }
  ret$improve <- NA
  if(!is.null(max.cand)) ret$improve[rightIndex %in% forward$rightIndex] <- forward$improve
  if(family == "gaussKern") {
    n <- length(y)
    # jump positions
    r <- ret$rightIndex
    ir <- r[-length(r)] # inner positions
    # correlate y and s, result is shorter than y
    ys <- convolve(y, s, conj = TRUE, type = "filter")
    # X' y is made up of left step up, constant centre, right step down
    ret$lXy <- c(c(rep(NA, kj - 1), ys)[ir], 0) # if jump has been found at kern$jump it corresponds to ys[1]
    ret$lcXy <- c(c(cumSum, rep(NA, kl - kj))[ir + kl - kj], NA)
    ret$rcXy <- c(rep(NA, kj), cumSum)[c(ir, n + kj)]
    ret$rXy <- c(c(cumSum, rep(NA, kl - kj))[ir + kl - kj] - c(rep(NA, kj), cumSum)[ir], 0) - ret$lXy # response to step from 1 to 0 is equal to 1 - step response
#     # check whether sparse matrices can be used
#     attr(ret, "Matrix") <- require(Matrix, quietly = TRUE)
  }
#   if(family == "gaussBound") {
#     # collapse bounds
#     ret$lower <- as.numeric(sapply(1:length(rightIndex), function(i) max(lower[c(1, rightIndex)[i]:rightIndex[i]])))
#     ret$upper <- as.numeric(sapply(1:length(rightIndex), function(i) min(upper[c(1, rightIndex)[i]:rightIndex[i]])))
#   }
  rownames(ret) <- 1:nrow(ret)
  class(ret) <- c("stepcand", class(ret))
  ret
}
