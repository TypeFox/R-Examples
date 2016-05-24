### maximal intersections from an interval-type Surv object
maximalint <- function(x, eps = 1e-7)
{
  if (is.null(nrow(x))) x <- matrix(x, nrow = 1)
  if (ncol(x) == 2) x[is.na(x[, 2]), 2] <- Inf
  else if (ncol(x) > 2) {
    ## right-censored
    x[x[, 3] == 0, 2] <- Inf
    ## exact
    x[x[, 3] == 1, 2] <- x[x[, 3] == 1, 1]
  }
  ind <- x[, 1] == x[, 2]
  ## break ties
  if (sum(!ind)) {
    l <- unique(x[!ind, 1])
    r <- r.minus <- unique(x[!ind, 2])
    r.minus[r %in% l] <- r[r %in% l] - eps
  }
  else l <- r <- r.minus <- NULL
  if (sum(ind)) {
    ## open left-endpoint
    x[ind, 1] <- x[ind, 2] - eps
    l <- c(l, unique(x[ind, 1]))
    r.minus <- c(r.minus, unique(x[ind, 2]))
    r <- c(r, unique(x[ind, 2]))
  }
  s <- rbind(cbind(l, l, 0), cbind(r.minus, r, 1))
  s <- s[order(s[, 1]), ]
  s <- cbind(s, rbind(s[-1, ], NA))
  ## maximal intersection left endpoint, right endpoint without and with ties 
  s <- s[s[, 3] == 0 & s[, 6] == 1, c(1, 4, 5)]
  if (is.null(nrow(s))) s <- matrix(s, nrow = 1)
  colnames(s) <- c("left", "right.untied", "right")
  ## maximal intersection overlap with censoring interval indicator matrix
  if (nrow(s) < 2) i <- matrix(1)
  else
    i <-
      t(apply(x[, 1:2], 1, function(x) 1 * (s[, 1] >= x[1] & s[, 2] <= x[2])))
  list(int = s, ind = i)
}
