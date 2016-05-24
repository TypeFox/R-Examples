CDVineSim <- function(N, family, par, par2 = rep(0, length(family)), type) {
  ## #############################################################################
  ##  Function that simulates copulae, gaussian, student's t, clayton or gumbel ##
  ## -------------------------------------------------------------------------- ##
  ##  INPUT:                                                                    ##
  ##    N                            Number of simulations                      ##
  ##    d                            Dimension of copula                        ##
  ##    family                       Array defining the ppc                     ##
  ##    par  		       Dependency parameters                      ##
  ##    par2			       Second copula parameter			  ##
  ##	  type			       Vine type (1=C-vine, 2=D-vine)		  ##
  ##  OUTPUT:                                                                   ##
  ##    U                            Simulated dxN copula                       ##
  ## -------------------------------------------------------------------------- ##
  ## #############################################################################
  
  if (type == "CVine") 
    type <- 1 else if (type == "DVine") 
    type <- 2
  if (type != 1 & type != 2) 
    stop("Vine model not implemented.")
  
  dd <- length(family)
  d <- (1 + sqrt(1 + 8 * dd))/2
  if (d != floor(d)) 
    stop("The length of the family vector is not correct.")
  
  if (length(par) < dd) 
    stop("Length of 'par' has to be d*(d-1)/2.")
  if (c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40) %in% family && length(par2) != 
    dd) 
    stop("Length of 'par2' has to be d*(d-1)/2 if t, BB1, BB6, BB7 or BB8 copulas are used.")
  if (length(family) != dd) 
    stop("Length of 'family' has to be d*(d-1)/2.")
  if (d < 2) 
    stop("Dimension has to be at least 2.")
  
  for (i in 1:(d * (d - 1)/2)) {
    if (!(family[i] %in% c(0, 1:10, 13, 14, 16:20, 23, 24, 26:30, 33, 34, 36:40))) 
      stop("Copula family not implemented.")
    # Parameterbereiche abfragen
    if ((family[i] == 1 || family[i] == 2) && abs(par[i]) >= 1) 
      stop("The parameter of the Gaussian and t-copula has to be in the interval (-1,1).")
    if (family[i] == 2 && par2[i] <= 2) 
      stop("The degrees of freedom parameter of the t-copula has to be larger than 2.")
    if ((family[i] == 3 || family[i] == 13) && par[i] <= 0) 
      stop("The parameter of the Clayton copula has to be positive.")
    if ((family[i] == 4 || family[i] == 14) && par[i] < 1) 
      stop("The parameter of the Gumbel copula has to be in the interval [1,oo).")
    if ((family[i] == 6 || family[i] == 16) && par[i] <= 1) 
      stop("The copula parameter of the Joe copula has to be in the interval (1,oo).")
    if (family[i] == 5 && par[i] == 0) 
      stop("The parameter of the Frank copula has to be unequal to 0.")
    if ((family[i] == 7 || family[i] == 17) && par[i] <= 0) 
      stop("The first parameter of the BB1 copula has to be positive.")
    if ((family[i] == 7 || family[i] == 17) && par2[i] < 1) 
      stop("The second parameter of the BB1 copula has to be in the interval [1,oo).")
    if ((family[i] == 8 || family[i] == 18) && par[i] <= 0) 
      stop("The first parameter of the BB6 copula has to be in the interval [1,oo).")
    if ((family[i] == 8 || family[i] == 18) && par2[i] < 1) 
      stop("The second parameter of the BB6 copula has to be in the interval [1,oo).")
    if ((family[i] == 9 || family[i] == 19) && par[i] < 1) 
      stop("The first parameter of the BB7 copula has to be in the interval [1,oo).")
    if ((family[i] == 9 || family[i] == 19) && par2[i] <= 0) 
      stop("The second parameter of the BB7 copula has to be positive.")
    if ((family[i] == 10 || family[i] == 20) && par[i] < 1) 
      stop("The first parameter of the BB8 copula has to be in the interval [1,oo).")
    if ((family[i] == 10 || family[i] == 20) && (par2[i] <= 0 || par2[i] > 1)) 
      stop("The second parameter of the BB8 copula has to be in the interval (0,1].")
    if ((family[i] == 23 || family[i] == 33) && par[i] >= 0) 
      stop("The parameter of the rotated Clayton copula has to be negative.")
    if ((family[i] == 24 || family[i] == 34) && par[i] > -1) 
      stop("The parameter of the rotated Gumbel copula has to be in the interval (-oo,-1].")
    if ((family[i] == 26 || family[i] == 36) && par[i] >= -1) 
      stop("The parameter of the rotated Joe copula has to be in the interval (-oo,-1).")
    if ((family[i] == 27 || family[i] == 37) && par[i] >= 0) 
      stop("The first parameter of the rotated BB1 copula has to be negative.")
    if ((family[i] == 27 || family[i] == 37) && par2[i] > -1) 
      stop("The second parameter of the rotated BB1 copula has to be in the interval (-oo,-1].")
    if ((family[i] == 28 || family[i] == 38) && par[i] >= 0) 
      stop("The first parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
    if ((family[i] == 28 || family[i] == 38) && par2[i] > -1) 
      stop("The second parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
    if ((family[i] == 29 || family[i] == 39) && par[i] > -1) 
      stop("The first parameter of the rotated BB7 copula has to be in the interval (-oo,-1].")
    if ((family[i] == 29 || family[i] == 39) && par2[i] >= 0) 
      stop("The second parameter of the rotated BB7 copula has to be negative.")
    if ((family[i] == 30 || family[i] == 40) && par[i] > -1) 
      stop("The first parameter of the rotated BB8 copula has to be in the interval (-oo,-1].")
    if ((family[i] == 30 || family[i] == 40) && (par2[i] >= 0 || par2[i] < (-1))) 
      stop("The second parameter of the rotated BB8 copula has to be in the interval [-1,0).")
  }
  
  if (type == 1) {
    if (!any(family %in% c(2, 7:10, 17:20, 27:30, 37:40))) 
      tmp <- .C("pcc", as.integer(N), as.integer(d), as.integer(family), as.integer(type), as.double(par), 
        as.double(rep(0, d * (d - 1)/2)), as.double(rep(0, N * d)), PACKAGE = "CDVine")[[7]] else tmp <- .C("pcc", as.integer(N), as.integer(d), as.integer(family), as.integer(type), as.double(par), 
      as.double(par2), as.double(rep(0, N * d)), PACKAGE = "CDVine")[[7]]
    U <- matrix(tmp, ncol = d)
    
  } else if (type == 2) {
    # D-vine
    
    Matrix <- matrix(rep(0, d * d), d, d)
    Copula.Params <- matrix(rep(0, d * d), d, d)
    Copula.Params2 <- matrix(rep(0, d * d), d, d)
    Copula.Types <- matrix(rep(0, d * d), d, d)
    
    for (i in 1:d) {
      Matrix[(d - i + 1), (d - i + 1)] <- i
    }
    
    k <- 1
    for (i in 1:(d - 1)) {
      for (j in 1:(d - i)) {
        Matrix[(d - i + 1), (d - j - i + 1)] <- j
        Copula.Types[(d - i + 1), (d - j - i + 1)] <- family[k]
        Copula.Params[(d - i + 1), (d - j - i + 1)] <- par[k]
        Copula.Params2[(d - i + 1), (d - j - i + 1)] <- par2[k]
        k <- k + 1
      }
    }
    
    Matrix[is.na(Matrix)] <- 0
    Copula.Types[is.na(Copula.Types)] <- 0
    Copula.Params[is.na(Copula.Params)] <- 0
    Copula.Params2[is.na(Copula.Params2)] <- 0
    
    MaxMat <- createMaxMat(Matrix)
    CondDistr <- neededCondDistr(Matrix)
    
    dvlist <- list(Matrix = Matrix, family = Copula.Types, par = Copula.Params, par2 = Copula.Params2, 
      MaxMat = MaxMat, CondDistr = CondDistr)
    
    o <- diag(dvlist$Matrix)
    dvlist <- normalizeDv(dvlist)
    
    matri <- as.vector(dvlist$Matrix)
    w1 <- as.vector(dvlist$family)
    th <- as.vector(dvlist$par)
    th2 <- as.vector(dvlist$par2)
    maxmat <- as.vector(dvlist$MaxMat)
    conindirect <- as.vector(dvlist$CondDistr$indirect)
    matri[is.na(matri)] <- 0
    w1[is.na(w1)] <- 0
    th[is.na(th)] <- 0
    th2[is.na(th2)] <- 0
    maxmat[is.na(maxmat)] <- 0
    conindirect[is.na(conindirect)] <- 0
    
    tmp <- rep(0, d * N)
    
    tmp <- .C("SimulateVine", as.integer(N), as.integer(d), as.integer(w1), as.integer(maxmat), as.integer(matri), 
      as.integer(conindirect), as.double(th), as.double(th2), as.double(tmp), PACKAGE = "CDVine")[[9]]
    
    out <- matrix(tmp, ncol = d)
    
    
    U <- out[, sort(o[length(o):1], index.return = TRUE)$ix]
    
  }
  return(U)
}

normalizeDv <- function(dvlist) {
  
  oldOrder <- diag(dvlist$Matrix)
  Matrix <- reorderDv(dvlist$Matrix)
  
  dvlist$Matrix[is.na(dvlist$Matrix)] <- 0
  dvlist$family[is.na(dvlist$family)] <- 0
  dvlist$par[is.na(dvlist$par)] <- 0
  dvlist$par2[is.na(dvlist$par2)] <- 0
  
  MaxMat <- createMaxMat(Matrix)
  CondDistr <- neededCondDistr(Matrix)
  
  return(list(Matrix = Matrix, family = dvlist$family, par = dvlist$par, par2 = dvlist$par2, MaxMat = MaxMat, 
    CondDistr = CondDistr))
}

reorderDv <- function(Matrix) {
  oldOrder <- diag(Matrix)
  
  O <- apply(t(1:nrow(Matrix)), 2, "==", Matrix)
  
  for (i in 1:nrow(Matrix)) {
    Matrix[O[, oldOrder[i]]] <- nrow(Matrix) - i + 1
  }
  
  return(Matrix)
}

createMaxMat <- function(Matrix) {
  
  if (dim(Matrix)[1] != dim(Matrix)[2]) 
    stop("Structure matrix has to be quadratic.")
  
  MaxMat <- reorderDv(Matrix)
  
  n <- nrow(MaxMat)
  
  for (j in 1:(n - 1)) {
    for (i in (n - 1):j) {
      MaxMat[i, j] <- max(MaxMat[i:(i + 1), j])
    }
  }
  
  tMaxMat <- MaxMat
  tMaxMat[is.na(tMaxMat)] <- 0
  
  oldSort <- diag(Matrix)
  oldSort <- oldSort[n:1]
  
  for (i in 1:n) {
    MaxMat[tMaxMat == i] <- oldSort[i]
  }
  
  return(MaxMat)
}

neededCondDistr <- function(Vine) {
  
  if (dim(Vine)[1] != dim(Vine)[2]) 
    stop("Structure matrix has to be quadratic.")
  
  Vine <- reorderDv(Vine)
  
  MaxMat <- createMaxMat(Vine)
  
  d <- nrow(Vine)
  
  M <- list()
  M$direct <- matrix(FALSE, d, d)
  M$indirect <- matrix(FALSE, d, d)
  
  M$direct[2:d, 1] <- TRUE
  
  for (i in 2:(d - 1)) {
    v <- d - i + 1
    
    bw <- as.matrix(MaxMat[i:d, 1:(i - 1)]) == v
    
    direct <- Vine[i:d, 1:(i - 1)] == v
    
    M$indirect[i:d, i] <- apply(as.matrix(bw & (!direct)), 1, any)
    
    M$direct[i:d, i] <- TRUE
    
    M$direct[i, i] <- any(as.matrix(bw)[1, ] & as.matrix(direct)[1, ])
  }
  
  return(M)
} 
