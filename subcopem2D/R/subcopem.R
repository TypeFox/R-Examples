subcopem <-
function(mat.xy, display = FALSE){
# 
# Input:  mat.xy = 2-column matrix with bivariate observations of (X,Y)
#                  ** repeated values are allowed
#                  ** slow for continuous sample size > 2000 (time > 2 min)
#                     => should better use <subcopemc> instead
#
# Output: Empirical 2-subcopula matrix, induced partitions,
#         standardized sample, and dependence measures
#              depMon = monotone standardized L1 distance [-1,1]
#               depL1 = standardized L1 distance [0,1]
#               depL2 = standardized L2 distance [0,1]
#         ** if display = TRUE dep measures and graphs are displayed
#
  n <- nrow(mat.xy)  # sample size
  X <- mat.xy[ , 1]  # X observed values
  Y <- mat.xy[ , 2]  # Y observed values
# Subcopula:
  R <- sort(unique(X))  # unique values of X (sorted)
  m1 <- length(R)  # number of unique values of X
  S <- sort(unique(Y))  # unique values of Y (sorted)
  m2 <- length(S)  # number of unique values of Y
  contar.aux <- function(r, s) sum((X <= r)*(Y <= s))  # counting function
  contar <- function(r, s) mapply(contar.aux, r, s)  # vectorized counting function
  subcopula <- matrix(0, nrow = (m1 + 1), ncol = (m2 + 1))  # creating subcopula matrix
  subcopula[2:(m1 + 1), 2:(m2 + 1)] <- (1/n)*outer(R, S, contar)  # subcopula calculation
# Partitions:
  p1 <- as.vector(table(sort(X))/n)  # observed proportion of unique values of X
  p2 <- as.vector(table(sort(Y))/n)  # observed proportion of unique values of Y
  q1 <- cumsum(p1)  # cumulative values of p1
  q2 <- cumsum(p2)  # cumulative values of p2
  partQ1 <- c(0, q1)  # first partition
  partQ2 <- c(0, q2)  # second partiion
# Standardized sample:
  f.aux <- function(x, A, q) q[which(x == A)]  # mapping function A -> q
  fR <- function(x) sapply(x, f.aux, A = R, q = q1)  # mapping function for X
  fS <- function(y) sapply(y, f.aux, A = S, q = q2)  # mapping function for Y
  muestra.std <- matrix(nrow = n, ncol = 2)  # creating matrix
  muestra.std[ , 1] <- fR(mat.xy[ , 1])  # calculating standardized sample
  muestra.std[ , 2] <- fS(mat.xy[ , 2])  # calculating standardized sample
# Dependence measures:
  M <- function(u, v) (u + v - abs(u - v))/2  # upper bound copula
  W <- function(u, v) (u + v - 1 + abs(u + v - 1))/2  # lower bound copula
  P <- function(u, v) u*v  # product copula (independence)
  subcopM <- outer(partQ1, partQ2, M)  # upper bound subcopula
  subcopW <- outer(partQ1, partQ2, W)  # lower bound subcopula
  subcopP <- outer(partQ1, partQ2, P)  # product subcopula (independence)
  dsgn <- function(A, B) sum(A - B)  # signed distance
  dL1 <- function(A, B) sum(abs(A - B))  # L1 distance
  dL2 <- function(A, B) sqrt(sum((A - B)^2))  # L2 distance
  Msgn <- dsgn(subcopM, subcopP)  # upper bound signed distance
  Wsgn <- dsgn(subcopP, subcopW)  # lower bound signed distance
  sgnBound <- max(Msgn, Wsgn)  # signed distance bound
  depMon <- dsgn(subcopula, subcopP)/sgnBound  # monotone standardized L1 distance
  ML1 <- dL1(subcopM, subcopP)  # upper bound L1 distance
  WL1 <- dL1(subcopP, subcopW)  # lower bound L1 distance
  L1bound <- max(ML1, WL1)  # L1 distance bound
  depL1 <- dL1(subcopula, subcopP)/L1bound  # standardized L1 distance
  ML2 <- dL2(subcopM, subcopP)  # upper bound L2 distance
  WL2 <- dL2(subcopP, subcopW)  # lower bound L2 distance
  L2bound <- max(ML2, WL2)  # L2 distance bound
  depL2 <- dL2(subcopula, subcopP)/L2bound  # L2 distance bound
# Output:
  SC <- list(depMon = depMon, depL1 = depL1, depL2 = depL2, 
             matrix = subcopula, part1 = partQ1, part2 = partQ2,
             sample.size = n, std.sample = muestra.std, sample = mat.xy)
  if (display == TRUE){
    message("monotone distance = ", round(depMon, 8))
    message("  std L1 distance = ", round(depL1, 8))
    message("  std L2 distance = ", round(depL2, 8))
    dev.new(); par(mfcol = c(2, 3))
    hist(mat.xy[ , 1], main = "histogram", xlab = "X")
    hist(mat.xy[ , 2], main = "histogram", xlab = "Y")
    plot(mat.xy, main = "sample", xlab = "X", ylab = "Y")
    plot(c(0, 1), c(0, 1), type = "n", main = "standardized sample", ylab = "",
         xlab = paste("monotone distance =", round(depMon, 3)))
    points(muestra.std)
    contour(partQ1, partQ2, subcopula, nlevels = 20, main = "subcopula", ylab = "",
            xlab = paste("std L2 distance =", round(depL2, 3)))
    image(partQ1, partQ2, subcopula, col = heat.colors(20), main = "subcopula", ylab = "",
          xlab = paste("std L1 distance =", round(depL1, 3)))   
  }
  return(SC) 
}
