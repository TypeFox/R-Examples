"stepbound" <-
function (y, bounds, ...) UseMethod("stepbound")

"stepbound.default" <-
function (y, bounds, x = 1:length(y), x0 = 2 * x[1] - x[2], max.cand = NULL, family = c("gauss", "gaussvar", "poisson", "binomial", "gaussKern"), param = NULL, weights = rep(1, length(y)), refit = y, jumpint = confband, confband = FALSE, ...)
{
  family <- match.arg(family)
  cand <- stepcand(y = y, x = x, x0 = x0, max.cand = max.cand, family = family, param = param, weights = weights, cand.radius = Inf)
  if(missing(bounds)) bounds <- stepR::bounds(y, family = family, param = param, ...)
  sb <- stepbound.stepcand(cand, bounds, refit) #, jumpint = jumpint, confband = confband

  if(confband) jumpint <- TRUE
  if(jumpint) {
    bounds.rev <- bounds
    bounds.rev$bounds$li <- length(y) - bounds$bounds$ri + 1
    bounds.rev$bounds$ri <- length(y) - bounds$bounds$li + 1
    bounds.rev$bounds <- bounds.rev$bounds[order(bounds.rev$bounds$li, bounds.rev$bounds$ri),]
    bounds.rev$start <- rep(NA, length(y))
    w <- which(diff(c(0,bounds.rev$bounds$li)) > 0)
    bounds.rev$start[bounds.rev$bounds$li[w]] <- w - 1 # R style to C style indices
    sb.rev <- stepbound.default(rev(y), bounds.rev, max.cand = max.cand, family = family, param = param, weights = weights, jumpint = FALSE, confband = FALSE, ...) # compute for reversed    data
#     bounds$bounds$ri <- length(y) - bounds$bounds$ri + 1
#     bounds$bounds$li <- length(y) - bounds$bounds$li + 1
#     names(bounds$bounds)[names(bounds$bounds) == "ri"] <- "RI"
#     names(bounds$bounds)[names(bounds$bounds) == "li"] <- "ri"
#     names(bounds$bounds)[names(bounds$bounds) == "RI"] <- "li"
#     bounds$bounds <- bounds$bounds[order(bounds$bounds$li, bounds$bounds$ri),]
#     bounds$start <- c(0, which(diff(bounds$bounds$li) > 0))
#     sb.rev <- stepbound(rev(y), bounds, max.cand = max.cand, family = family, param = param, weights = weights, jumpint = FALSE, confband = FALSE, ...) # compute for reversed data
    sb$rightIndexLeftBound <- c(pmax(sb$rightIndexLeftBound[-nrow(sb)], rev(length(y) - sb.rev$rightIndexRightBound[-nrow(sb.rev)])), length(y)) # combining bounds from original and reversed data
    sb$rightIndexRightBound <- c(pmin(sb$rightIndexRightBound[-nrow(sb)], rev(length(y) - sb.rev$rightIndexLeftBound[-nrow(sb.rev)])), length(y))
    sb$leftIndexRightBound <- c(1, sb$rightIndexRightBound[-nrow(sb)]+1)
    sb$leftIndexLeftBound <- c(1, sb$rightIndexLeftBound[-nrow(sb)]+1)
    sb$rightEndLeftBound <- x[sb$rightIndexLeftBound]
    sb$rightEndRightBound <- x[sb$rightIndexRightBound]
    sb$leftEndLeftBound <- x[sb$leftIndexLeftBound]
    sb$leftEndRightBound <- x[sb$leftIndexRightBound]
  }
  
  if(confband) {
    band <- .Call(.confBand, as.integer(sb$rightIndexLeftBound), as.integer(c(0, sb$rightIndexRightBound[-nrow(sb)])), as.integer(bounds$start), as.integer(bounds$bounds$ri - 1), as.numeric(bounds$bounds$lower), as.numeric(bounds$bounds$upper))
    band <- cbind(x = x, as.data.frame(band))
    attr(band, "x0") <- x0
    class(band) <- c("confband", class(band))
    attr(sb, "confband") <- band
  }
  
  sb
}

"stepbound.stepcand" <-
function (y, bounds, refit = TRUE, ...)
{
  # ensure start matches indices of y
  if(nrow(y) != length(bounds$start)) {
    stop("length(bounds$start) ", length(bounds$start), " and nrow(y) ", nrow(y), " don't match")
#     if(max(y$rightIndex) != length(bounds$start))
#     stop("bounds$start and y$rightIndex don't match")
#     # select the bounds corresponding to the candidates
#     start <- bounds$start[y$rightIndex]
#     end <- c(bounds$start[y$rightIndex[-nrow(y)] + 1] - 1, bounds$start[y$rightIndex[nrow(y)]])
#     print(y$rightIndex)
#     print(start)
#     print(end)
#     bounds$start <- c(1, cumsum(end - start + 1)[-length(start)]) - 1
#     bounds$bounds <- bounds$bounds[unlist(apply(rbind(start, end), 2, function(i) i[1]:i[2])),]
#     print(bounds)
  }
  # check feasibility
  if(!bounds$feasible) stop("Bounds not feasible!")
  algo <- switch(attr(y, "family"),
    gaussKern = "gauss",
    attr(y, "family")
  )
  b <- switch(algo,
    gauss = {
      .Call(.boundedGauss, y$cumSum, y$cumSumSq, y$cumSumWe, as.integer(bounds$start), as.integer(bounds$bounds$ri - 1), as.numeric(bounds$bounds$lower), as.numeric(bounds$bounds$upper))
    },
    gaussvar = {
      # assume that both bounds are 0 where y is constant 0
      .Call(.boundedGaussVar, y$cumSumSq, y$cumSumWe, as.integer(bounds$start), as.integer(bounds$bounds$ri - 1), as.numeric(bounds$bounds$lower), as.numeric(bounds$bounds$upper))
    },
    poisson = {
      if(any(bounds$bounds$upper < 0)) stop("negative upper bounds not permitted!")
      .Call(.boundedPoisson, y$cumSum, y$cumSumWe, as.integer(bounds$start), as.integer(bounds$bounds$ri - 1), as.numeric(bounds$bounds$lower), as.numeric(bounds$bounds$upper))
    },
    binomial = {
      if(any(bounds$bounds$upper < 0)) stop("negative upper bounds not permitted!")
      if(any(bounds$bounds$lower > 1)) stop("lower bounds > 1 not permitted!")
      .Call(.boundedBinom, as.integer(attr(y, "param")), y$cumSum, y$cumSumWe, as.integer(bounds$start), as.integer(bounds$bounds$ri - 1), as.numeric(bounds$bounds$lower), as.numeric(bounds$bounds$upper))
    },
    stop("unknown family")
  )
  ret <- stepfit(cost = attr(b, "cost"), family = attr(y, "family"), value = b$value, param = attr(y, "param"), 
    leftEnd = y$leftEnd[c(0, b$rightEnd[-length(b$rightEnd)]) + 1], rightEnd = y$rightEnd[b$rightEnd], x0 = attr(y, "x0"),
    leftIndex = y$leftIndex[c(0, b$rightEnd[-length(b$rightEnd)]) + 1], rightIndex = y$rightIndex[b$rightEnd])
  ret$rightIndexLeftBound <- y$rightIndex[b$endLeftBound]
  ret$rightIndexRightBound <- y$rightIndex[b$endRightBound]
  ret$rightEndLeftBound <- y$rightEnd[b$endLeftBound]
  ret$rightEndRightBound <- y$rightEnd[b$endRightBound]
  if(attr(y, "family") != "gaussvar") {
    ret$cumSum <- y$cumSum[b$rightEnd]
  }
  ret$cumSumWe <- y$cumSumWe[b$rightEnd]
  if(attr(y, "family") %in% c("gauss", "gaussvar")) {
    ret$cumSumSq <- y$cumSumSq[b$rightEnd]
  }
  if(attr(y, "family") == "gaussKern") {
    ret$cumSumSq <- y$cumSumSq[b$rightEnd]
    ret$lXy <- y$lXy[b$rightEnd]
    ret$lcXy <- y$lcXy[b$rightEnd]
    ret$rcXy <- y$rcXy[b$rightEnd]
    ret$rXy <- y$rXy[b$rightEnd]
    if(!identical(refit, FALSE) & nrow(ret) > 1) ret <- ret[refit = refit]
  }

  ret
}
