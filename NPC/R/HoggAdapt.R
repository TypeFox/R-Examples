HoggAdapt <- function (y) {
  ## Function for adaptive choice of test based on skewness and tailweight
  ## Best used only for continuous distributions that differ in location
  ## only (not in scale)
  N <- length(y)
  Zinc <- sort(y, decreasing=FALSE) ## order statistics
  Zdec <- sort(y, decreasing=TRUE) ## order statistics
  M.5 <- mean(y, trim=.25)
  ## tails
  int.05 <- trunc(N / 20)
  decim.05 <- (N %% 20) / 20
  L.05 <- (sum(Zinc[1:int.05]) + decim.05*Zinc[int.05 + 1])/(N/20)
  U.05 <- (sum(Zdec[1:int.05]) + decim.05*Zdec[int.05 + 1])/(N/20)
  int.5 <- trunc(N / 2)
  decim.5 <- (N %% 2) / 2
  L.5 <- (sum(Zinc[1:int.5]) + decim.5*Zinc[int.5 + 1])/(N/2)
  U.5 <- (sum(Zdec[1:int.5]) + decim.5*Zdec[int.5 + 1])/(N/2)
  skew <- (U.05 - M.5)/(M.5 - L.05)
  tailwt <- (U.05 - L.05)/(U.5 - L.5)
  dist <- NULL
  if (tailwt > 7) {
    dist <- 'Very Heavy-Tailed Distribution'
    fun <- function (y, tr, tl) {
      ## Very heavy tails -> median test
      Ri <- rank(y)
      Ai <- as.integer(Ri > (N + 1)/2)
      T <- mean.default(Ai[tr==tl]) - mean.default(Ai[tr!=tl])
      return(T)
    }
  }
  if (skew > 2 && tailwt <= 7) {
    dist <- 'Right-Skewed Distribution'
    fun <- function (y, tr, tl) {
      ## Right skewed -> drop top half of distribution
      Ri <- rank(y)
      Ai <- ifelse(Ri <= floor((N + 1)/2), Ri - floor((N + 1)/2), 0)
      T <- mean.default(Ai[tr==tl]) - mean.default(Ai[tr!=tl])
      return(T)
    }
  }
  if (skew < 1/2 && tailwt <= 7) {
    dist <- 'Left-Skewed Distribution'
    fun <- function (y, tr, tl) {
      ## Left skewed -> drop bottom half of distribution
      Ri <- rank(y)
      Ai <- ifelse(Ri >= floor((N + 1)/2), Ri - floor((N + 1)/2), 0)
      T <- mean.default(Ai[tr==tl]) - mean.default(Ai[tr!=tl])
      return(T)
    }
  }
  if (skew >= 1/2 && skew <= 2 && tailwt <= 2) {
    dist <- 'Light-Tailed Symmetric Distribution'
    fun <- function (y, tr, tl) {
      ## Light-tailed symmetric -> drop middle half of distribution
      Ri <- rank(y)
      Ai <- 0
      Ai <- ifelse(Ri <= floor((N + 1)/4),
                   Ri - floor((N + 1)/4) - 1/2, Ai)
      Ai <- ifelse(Ri >= N - floor((N + 1)/4) + 1,
                   Ri - N + floor((N + 1)/4) - 1/2, Ai)
      T <- mean.default(Ai[tr==tl]) - mean.default(Ai[tr!=tl])
      return(T)
    }
  }
  if (is.null(dist)) {
    dist <- 'Heavy-Tailed Distribution'
    fun <- function (y, tr, tl) {
      ## Heavy tails -> regular rank test
      Ri <- rank(y)
      Ai <- Ri
      T <- mean.default(Ai[tr==tl]) - mean.default(Ai[tr!=tl])
      return(T)
    }
  }
  return(fun)
}
