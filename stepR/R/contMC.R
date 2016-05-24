# simulate continuous time Markov chain, along the lines of
# Antonius M. J. VanDongen (1996) "A New Algorithm for Idealizing Single Ion Channel Data Containing Multiple Unknown Conductance Levels", Biophysical Journal 70, 1303-1315.
"contMC" <-
function(n, values, rates, start = 1, sampling = 1, family = c("gauss", "gaussKern"), param = NULL)
{
  family <- match.arg(family)
  # rates[i,j] specifies rate of leaving state i and going to state j
  stopifnot(rep(length(values), 2) == dim(rates))
  
  # rate of leaving state i
  diag(rates) <- 0
  leave <- rowSums(rates)
  stopifnot(all(leave > 0))
  
  # length of sequence to simulate
  N <- n
  L <- n / sampling
  if(family == "gauss") {
    stopifnot(is.numeric(param))
    stopifnot(param >= 0)
  }
  if(family == "gaussKern") {
    stopifnot("dfilter" %in% class(param$df))
    stopifnot(is.numeric(param$sd))
    stopifnot(param$sd >= 0)
    # over sampling
    stopifnot(is.numeric(param$over))
    param$over <- as.integer(param$over)
    stopifnot(param$over > 0)
    N <- param$over * n + length(param$df$kern) - 1
    L <- L + ( length(param$df$kern) - 1 ) / param$over / sampling
  }
  
  # simulate events
  s <- 0
  len <- list()
  state <- list(start)
  i <- 1
  while(s < L) {
    len[[i]] <- rexp(1, leave[state[[i]]])
    s <- s + len[[i]]
    i <- i + 1
    state[[i]] <- sample(length(values), 1, prob = rates[state[[i - 1]],])
  }
  len <- unlist(len)
  state <- unlist(state)[-i]
  re <- cumsum(len)
  
  cont <- stepblock(values[state], leftEnd = c(0, re[-length(re)]), rightEnd = re)
  cont$state <- state
  
  # sample
  sr <- if(family == "gaussKern") sampling * param$over else sampling
  ri <- floor(re * sr)
  # remove events shorter than sampling interval
  keep <- which(c(diff(ri) > 0, TRUE))
  ri <- ri[keep]
  ri[length(ri)] <- N
  li <- c(0, ri[-length(ri)]) + 1
  state <- state[keep]
  
  discr <- stepfit(NA, family, values[state], param, li / sr, ri / sr, 0, li, ri)
  discr$states <- state
  
  # add noise
  y <- switch(family,
    gauss = fitted(discr) + rnorm(N, 0, param),
    gaussKern = {
      attr(discr, "family") <- "gauss"
      f <- fitted(discr)
      attr(discr, "family") <- "gaussKern"
      f + rnorm(N, 0, param$sd / sqrt(sum(param$df$kern^2)))
    },
    stop("family not implemented yet")
  )
  
  if(family == "gaussKern") {
    # filter
    y <- filter(y, rev(param$df$kern), "convolution", 1) # first length(param$df$kern) - 1 values are NA
    # down-sample and align
    y <- y[param$over * (0:(n-1)) + length(param$df$kern)]
    cont$rightEnd <- cont$rightEnd - param$df$jump / sr + 0.5 / sampling
    cont$leftEnd <- cont$leftEnd - param$df$jump / sr + 0.5 / sampling
    ri <- pmax(0, floor(0.5 + ( ri - param$df$jump ) / param$over))
    keep <- which(c(diff(ri) > 0, TRUE) & c(TRUE, ri[-length(ri)] < n))
    ri <- ri[keep]
    ri[length(ri)] <- n
    li <- c(0, ri[-length(ri)]) + 1
    state <- state[keep]
    discr <- stepfit(NA, family, values[state], param, li / sampling, ri / sampling, 0, li, ri)
    discr$states <- state
  }
  
  list(cont = cont, discr = discr, data = data.frame(x = (1:n) / sampling, y = y))
}
