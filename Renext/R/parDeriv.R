##******************************************************************
## Author: Yves Deville
## 
## Exact derivation of the log-likelihood and the
## distribution function of the folowing models
##
##   1) exponential [rate]                  NOT IMPLEMENTED
##   2) weibull     [shape, scale]
##   3) GPD         [scale, shape]
##   4) mixture of two exponentials         NOT IMPLEMENTED
## 
## derlogdens, derlogdens2 : derivatives of log-density log f(s) 
## derSurv, der2Surv:        derivatives of the survival S(x)   
## 
##****************************************************************** 

parDeriv <- function(par, x, distname, sum = TRUE) {
  
  n <- length(x)
 
  if (distname == "weibull") {
    
    parnames <- c("shape", "scale")
    alpha <- par["shape"]
    beta <- par["scale"]
    
    x1 <- x / beta
    U <- log(x1)
    H <- x1^alpha
    S <- exp(-H)

    derH <- matrix(rep(0, n * 2L), nrow = length(x), ncol = 2L)
    derU <- matrix(rep(0, n * 2L), nrow = length(x), ncol = 2L)  
    colnames(derH) <- colnames(derU) <- parnames

    UH <- U * H
    derH[ ,"shape"] <- UH
    derH[ , "scale"] <- -alpha * H / beta
    derU[ , "scale"] <- -1 / beta

    derlogdens <- matrix(rep(0, length(x)*2L), nrow = length(x), ncol = 2L)
    colnames(derlogdens) <- parnames
    
    derlogdens[ , "shape"] <- 1 / alpha + U - UH
    derlogdens[ , "scale"] <- alpha * (-1 + H) / beta
    
    der2H <- array(0, dim = c(n, 2L, 2L),
                   dimnames <- list(1L:n, parnames, parnames))
    der2logdens <- array(0, dim = c(n, 2L, 2L),
                         dimnames <- list(1L:n, parnames, parnames))
    
    der2H[ , "shape", "shape"] <- U * UH
    der2H[ , "shape", "scale"] <- -H * (1 + alpha * U) / beta
    der2H[ , "scale", "shape"] <-  der2H[ , "shape", "scale"] 
    der2H[ , "scale", "scale"] <- alpha * (alpha + 1) * H / beta / beta

    der2logdens[ , "shape", "shape"] <- -1 / alpha / alpha - U * UH
    der2logdens[ , "shape", "scale"] <- (-1 + H + alpha * UH) / beta
    der2logdens[ , "scale", "shape"] <- der2logdens[ , "shape", "scale"]
    der2logdens[ , "scale", "scale"] <- (1 - H - alpha * H) * alpha / beta / beta
      
    derF <- sweep(x = derH, MARGIN = 1, STATS = S, FUN = "*")
    der2F <- sweep(x = der2H, MARGIN = 1, STATS = S, FUN = "*")
    Add <- t(apply(derH, 1, function(x) as.vector(x %o% x)))
    Add <- sweep(x = Add, MARGIN = 1, STATS = S, FUN = "*")
    Add <- array(Add, dim = c(n, 2L, 2L))
    der2F <- der2F - Add
    
  } else if (distname == "gpd") {
    
    parnames <- c("scale", "shape")

    xi <- par["shape"]
    xi1 <- xi + 1
    sigma <- par["scale"]
    z <- x / sigma
    
    A <- log(1 + xi*z)
    B <- z / (1 + xi*z)
    B2 <- B * B
    S <- exp(-A / xi)
    
    derA <- matrix(rep(0, n*2L), nrow = length(x), ncol = 2L)
    derB <- matrix(rep(0, n*2L), nrow = length(x), ncol = 2L)  
    colnames(derA) <- colnames(derB) <- parnames
    
    derA[ , "scale"] <- -xi * B / sigma
    derA[ ,"shape"] <- B
    
    derB[ , "scale"] <- -(B - xi * B2) / sigma
    derB[ , "shape"] <- -B2
    
    derlogdens <- matrix(rep(0, length(x)*2), nrow = length(x), ncol = 2)
    colnames(derlogdens) <- parnames
    
    derlogdens[ , "scale"] <- -(1 - (xi+1) * B) / sigma
    derlogdens[ , "shape"] <- A  / xi / xi - (xi + 1) * B / xi
    
    der2logdens <- array(0, dim = c(n, 2L, 2L),
                         dimnames <- list(1L:n, parnames, parnames))
    
    der2logdens[ , "scale", "scale"] <- (1 - 2 * xi1 * B + xi * xi1 * B2) / sigma / sigma
    der2logdens[ , "scale", "shape"] <- (B - xi1 * B2) / sigma
    der2logdens[ , "shape", "scale"] <- der2logdens[ , "scale", "shape"]
    der2logdens[ , "shape", "shape"] <- -2 * A / xi / xi / xi + 2 * B / xi / xi + xi1 * B2 / xi

    derH <- matrix(rep(0, length(x)*2), nrow = length(x), ncol = 2L)
    colnames(derH) <- parnames
  
    derH[ , "scale"] <- -B / sigma
    derH[ , "shape"] <- -A / xi / xi + B / xi
    
    der2H <- array(0, dim = c(n, 2L, 2L),
                   dimnames <- list(1L:n, parnames, parnames))
    
    der2H[ , "scale", "scale"] <- (2 * B - xi * B2) / sigma / sigma
    der2H[ , "scale", "shape"] <- B2/sigma
    der2H[ , "shape", "scale"] <- der2H[ , "scale", "shape"]
    der2H[ , "shape", "shape"] <- 2 * A / xi / xi / xi - 2 * B / xi / xi - B2 / xi

    derF <- sweep(x = derH, MARGIN = 1, STATS = S, FUN = "*")
    der2F <- sweep(x = der2H, MARGIN = 1, STATS = S, FUN = "*")
    Add <- t(apply(derH, 1, function(x) as.vector(x %o% x)))
    Add <- array(Add, dim = c(n, 2L, 2L))
    Add <- sweep(x = Add, MARGIN = 1, STATS = S, FUN = "*")
    der2F <- der2F - Add
   
  } else if (distname == "lomax") {
      
      parnames <- c("shape", "scale")

      alpha <- par["shape"]
      beta <- par["scale"]
      z <- 1 / (1 + x / beta) 
      S <- z^alpha
      
      derlogdens <- matrix(rep(0, length(x) * 2), nrow = length(x), ncol = 2)
      colnames(derlogdens) <- parnames
      
      derlogdens[ , "shape"] <- 1 / alpha + log(z)
      derlogdens[ , "scale"] <- ( alpha - (alpha + 1) * z) / beta
      
      der2logdens <- array(0, dim = c(n, 2L, 2L),
                           dimnames <- list(1L:n, parnames, parnames))
      
      der2logdens[ , "shape", "shape"] <- - 1 / alpha / alpha
      der2logdens[ , "shape", "scale"] <- (1 - z) / beta
      der2logdens[ , "scale", "shape"] <- der2logdens[ , "shape", "scale"]
      der2logdens[ , "scale", "scale"] <- (1 - (alpha + 1) * (1 - z^2)) / beta / beta
      
      derH <- matrix(rep(0, length(x) * 2), nrow = length(x), ncol = 2L)
      colnames(derH) <- parnames
      
      derH[ , "shape"] <- -log(z)
      derH[ , "scale"] <- - alpha * ( 1 - z) / beta
    
      der2H <- array(0, dim = c(n, 2L, 2L),
                   dimnames <- list(1L:n, parnames, parnames))
    
      der2H[ , "shape", "shape"] <- 0
      der2H[ , "shape", "scale"] <- -(1 - z) / beta
      der2H[ , "scale", "shape"] <- der2H[ , "shape", "scale"]
      der2H[ , "scale", "scale"] <- alpha * (1 - z^2) / beta / beta
      
      derF <- sweep(x = derH, MARGIN = 1, STATS = S, FUN = "*")
      der2F <- sweep(x = der2H, MARGIN = 1, STATS = S, FUN = "*")
      Add <- t(apply(derH, 1, function(x) as.vector(x %o% x)))
      Add <- array(Add, dim = c(n, 2L, 2L))
      Add <- sweep(x = Add, MARGIN = 1, STATS = S, FUN = "*")
      der2F <- der2F - Add
     

  } else if (distname == "maxlo") {
      
      parnames <- c("shape", "scale")

      alpha <- par["shape"]
      beta <- par["scale"]
      z <- 1 / (1 - x / beta) 
      S <- z^(-alpha)
      
      derlogdens <- matrix(rep(0, length(x) * 2), nrow = length(x), ncol = 2)
      colnames(derlogdens) <- parnames
      
      derlogdens[ , "shape"] <- 1 / alpha - log(z)
      derlogdens[ , "scale"] <- (-alpha + (alpha - 1) * z) / beta
      
      der2logdens <- array(0, dim = c(n, 2L, 2L),
                           dimnames <- list(1L:n, parnames, parnames))
      
      der2logdens[ , "shape", "shape"] <- - 1 / alpha / alpha
      der2logdens[ , "shape", "scale"] <- (z - 1) / beta
      der2logdens[ , "scale", "shape"] <- der2logdens[ , "shape", "scale"]
      der2logdens[ , "scale", "scale"] <- (1 - (alpha - 1) * (z^2 - 1)) / beta / beta
      
      derH <- matrix(rep(0, length(x) * 2), nrow = length(x), ncol = 2L)
      colnames(derH) <- parnames
      
      derH[ , "shape"] <- log(z)
      derH[ , "scale"] <- - alpha * (z - 1) / beta
    
      der2H <- array(0, dim = c(n, 2L, 2L),
                   dimnames <- list(1L:n, parnames, parnames))
    
      der2H[ , "shape", "shape"] <- 0
      der2H[ , "shape", "scale"] <- -(z - 1) / beta
      der2H[ , "scale", "shape"] <- der2H[ , "shape", "scale"]
      der2H[ , "scale", "scale"] <- alpha * (z^2 - 1) / beta / beta
      
      derF <- sweep(x = derH, MARGIN = 1, STATS = S, FUN = "*")
      der2F <- sweep(x = der2H, MARGIN = 1, STATS = S, FUN = "*")
      Add <- t(apply(derH, 1, function(x) as.vector(x %o% x)))
      Add <- array(Add, dim = c(n, 2L, 2L))
      Add <- sweep(x = Add, MARGIN = 1, STATS = S, FUN = "*")
      der2F <- der2F - Add
     

  }
  
  if (sum) {
    derlogdens <- apply(derlogdens, 2, sum)
    der2logdens <- apply(der2logdens, 2:3, sum)
    derF <- apply(derF, 2L, sum)
    der2F <- apply(der2F, 2L:3L, sum)
  }
  
  list(## der2H = der2H,
       derLogdens = derlogdens,
       der2Logdens = der2logdens,
       derSurv = -derF,
       der2Surv = -der2F)
  
  
}
