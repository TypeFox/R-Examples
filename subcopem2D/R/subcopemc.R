subcopemc <-
function(mat.xy, m = nrow(mat.xy), display = FALSE){
# 
# Input:  mat.xy = 2-column matrix with n bivariate observations of (X,Y)
#                  ** repeated values are NOT allowed
#                  ** slow if m > 2000 (time > 2 min)
#                  ** if display = TRUE dependence measures and graphs are displayed
#
# Output: Empirical 2-subcopula matrix of order m in {1,...,n}, 
#         induced partitions, standardized sample, and dependence measures
#              depMon = monotone standardized L1 distance [-1,1]
#               depL1 = standardized L1 distance [0,1]
#               depL2 = standardized L2 distance [0,1]
#         ** If m = n (default) it is the usual empirical subcopula
#
# Checking for appropriate order m
  n <- nrow(mat.xy)  # sample size
  if (!(m %in% (1:n))){
    error.msg <- paste("Order m must be an integer value in 1:n, n =", n)
    stop(error.msg)
  } else{
# Checking for repeated values:
  mensaje <- character(0)
  X <- mat.xy[ , 1]
  Y <- mat.xy[ , 2]
  if ((length(unique(X)) < n) | (length(unique(Y)) < n)){
    mensaje <- "Presence of repeated values, jittering has been applied"
    ind.X <- which(duplicated(X))
    ind.Y <- which(duplicated(Y))
    if (length(ind.X) > 0) X[ind.X] <- jitter(X[ind.X])
    if (length(ind.Y) > 0) Y[ind.Y] <- jitter(Y[ind.Y])
    mat.xy <- cbind(X, Y)
  }  
# Subcopula of order m:
  r.xy <- apply(mat.xy, 2, rank)  # bivariate ranks
  kparticion <- round(seq(0, n, length = (m + 1)), 0)  # integer partition of order m
  r.ordx <- r.xy[order(r.xy[ , 1]) , ]  # reordering bivariare ranks on first variable
  contar.aux <- function(kx, ky) sum(r.ordx[1:kx, 2] <= ky)  # counting function
  contar <- function(kx, ky) mapply(contar.aux, kx, ky)  # vectorized counting function
  subcopula <- matrix(0, nrow = (m + 1), ncol = (m + 1))  # creating subcopula matrix
  # subcopula calculation:
  subcopula[2:(m + 1), 2:(m + 1)] <- (1/n)*outer(kparticion[2:(m + 1)], kparticion[2:(m + 1)], contar)
# Induced partition:
  particion <- kparticion/n
# Standardized sample:
  muestra.std <- r.xy/n
# Dependence measures:
  M <- function(u, v) (u + v - abs(u - v))/2  # upper bound copula
  W <- function(u, v) (u + v - 1 + abs(u + v - 1))/2  # lower bound copula
  P <- function(u, v) u*v  # product copula (independence)
  subcopM <- outer(particion, particion, M)  # upper bound subcopula
  subcopW <- outer(particion, particion, W)  # lower bound subcopula
  subcopP <- outer(particion, particion, P)  # product subcopula (independence)
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
             matrix = subcopula, part1 = particion, part2 = particion,
             sample.size = n, order = m, std.sample = muestra.std, sample = mat.xy)
  if (length(mensaje) > 0) warning(mensaje)
  if (display == TRUE){
    message("sample size = ", n, "  order = ", m)
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
    contour(particion, particion, subcopula, nlevels = 20, main = "subcopula", ylab = "",
            xlab = paste("std L2 distance =", round(depL2, 3)))
    image(particion, particion, subcopula, col = heat.colors(20), main = "subcopula", ylab = "",
          xlab = paste("std L1 distance =", round(depL1, 3)))
  }
  return(SC)
  }
}
