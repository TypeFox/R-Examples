#### function for generating hypotheses matrices
# nf: number of factors involved
# fl: vector of factorlevels, i.e. (a, b, c, ....) of size nf

#--------------------- crossed designs ------------------------------
HC <- function(fl, perm_names, names){
  nf <- length(fl)
  # centering matrix
  P <- function(x){
    P <- diag(x) - matrix(1 / x, ncol = x, nrow = x)
    return(P)
  }
  # scaled one-vectors
  One <- function(x){
    I <- t(1 / x * matrix(1, ncol = 1, nrow = x))
    return(I)
  }
  # number of hypotheses
  tmp <- 0
  for (i in 1:nf) {
    tmp <- c(tmp, choose(nf, i))
    nh <- sum(tmp)
  }
  # calculate the permutation of the names
  Z <- 0:(nf - 1)
  position <- rep(0, nh)
  for (i in 1:nh) {
    position[i] <- position[i] + sum(2 ^ Z[which(perm_names[i, ] == 1)])
  }
  Vektor <- rep(NA, nh)
  for (j in 1:nh) {
    Vektor[position[j]] <- j
  }
  fac_names <- names[Vektor]
  # function for calculating the kronecker product of several matrices
  kp <- function(A) {
    kp <- A[[1]]
    for (i in 2: length(A)) {
      kp <- kp %x% A[[i]]
    }
    return(kp)
  }
  # scaled One-vectors for all factors
  A <- list()
  for (i in 1:nf) {
    A[[i]] <- One(fl[i])
  }
  A_alt <- A
  # P-matrices for all factors
  B <- list()
  for (i in 1:nf) {
    B[[i]] <- P(fl[i])
  }
  # calculate all combinations of elements of A and B
  # (except all from A and all from B)
  n1 <- c(0, 1)
  liste <- vector("list", nf)
  for (i in 1:nf) {
    liste[[i]] <- n1
  }
  G <- expand.grid(liste)
  G <- G[2:(dim(G)[[1]] - 1), ]        # 1 means A
  # list of all combinations of A and B
  C <- list()
  for (i in 1:dim(G)[1]) {
    index <- which(G[i, ] == 1)
    for (j in 1:length(index)) {
      A[[index[j]]] <- B[[index[j]]]
    }
    C[[i]] <- A
    A <- A_alt
  }
  # calculation of the hypotheses matrices (except for the nf-fold interaction)
  hypo <- vector("list", nh)
  for (i in 1:(nh - 1)) {
    C_tmp <- C[[i]]
    hypo[[i]] <- kp(C_tmp)
  }
  # nf-fold interaction
  hypo[[nh]] <- kp(B)
  return(list(hypo, fac_names))
}

# ----------------------- nested designs -------------------------------------
HN <- function(fl){
  nf <- length(fl)
  # centering matrix
  P <- function(x){
    P <- diag(x) - matrix(1 / x, ncol = x, nrow = x)
    return(P)
  }
  # scaled one-vectors
  One <- function(x){
    I <- t(1 / x * matrix(1, ncol = 1, nrow = x))
    return(I)
  }
  # function for calculating the kronecker product of several matrices
  kp <- function(A) {
    kp <- A[[1]]
    for (i in 2: length(A)) {
      kp <- kp %x% A[[i]]
    }
    return(kp)
  }
  hypo <- vector("list", nf)
  if (nf == 2) {
    fl[2] <- fl[2] / fl[1]
    hypo[[1]] <- P(fl[1]) %x% One(fl[2])
    hypo[[2]] <- diag(fl[1]) %x% P(fl[2])
  } else if (nf == 3) {
    fl[3] <- fl[3] / fl[2]
    fl[2] <- fl[2] / fl[1]
    hypo[[1]] <- P(fl[1]) %x% One(fl[2]) %x% One(fl[3])
    hypo[[2]] <- diag(fl[1]) %x% P(fl[2]) %x% diag(fl[3])
    hypo[[3]] <- diag(fl[1]) %x% diag(fl[2]) %x% P(fl[3])
  } else {
    print("Error")
  }
  return(hypo)
}
