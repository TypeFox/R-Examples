Bell <- function (n) {                                                           
  coeff <- vector("list", n)
  derivOrder <- vector("list", n)
  nways <- function(combi){
    ncol(setparts(combi))
  }
  for (i in 1 : n) {                                                             
    partOfn <- restrictedparts(n, i, include.zero = FALSE, 
      decreasing = FALSE)   
    k <- apply(partOfn, 2, nways) 
    if (i == 1) k <- 1 
    derivOrder[[i]] <- partOfn
    coeff[[i]] <- k
  }  
  out <- c(derivOrder, coeff)
  return(out)
}                       
                                                                              
derivErfc <- function(n, wrt1 = "a", wrt2 = "u", 
            expr2 = expression(u^2), aVal = 1, bVal = 2, 
            xVal = 1, varU = TRUE, lower = TRUE) {
  a <- aVal
  b <- bVal
  x <- xVal  
  u <- sqrt(b)/x + x * sqrt(a)           
  is.whole <- function(a) { 
    (is.numeric(a) && floor(a)==a) ||
    (is.complex(a) && floor(Re(a)) == Re(a) && 
      floor(Im(a)) == Im(a))
  }
  Delta <- function(n, expr = expression(- x * a^3), 
    wrt = "a", val = 1, xVal = 2) {
    derivatives <- function(order) { 
      d <- Deriv(expr = expr, wrt = wrt, order = order)  
      return(d)                               
    }   
  funCoeff <- function(n) {
    outCoeff <- vector('list', n)
    for (i in 1 : n) {
      k <- (n + i)
      outCoeff[[i]] <- Bell(n)[[n + i]]
    }
    out <- outCoeff 
    return(out)
  }
  funDerivOrder <- function(n) {
    outDerivOrder <- vector('list', n)
    nrowVec <- numeric(n)
    outCoeff <- vector('list', n)
    for (i in 1 : n) {
      outDerivOrder[[i]] <- Bell(n)[[i]]
      nrowVec[i] <- nrow(outDerivOrder[[i]]) 
    }
    out <- outDerivOrder
    return(out)
  }
  name <- "u"
  assign(name, wrt)
  wrt <- get(name)  
  MH <- numeric(n)
  assign(wrt, val)
  x <- xVal
  orders <- funDerivOrder(n)
  coeffs <- funCoeff(n) 
  for(i in 1:n) {
    nr <- nrow(orders[[i]])
    nc <- ncol(orders[[i]])
    derivs <- numeric(nr * nc)
    for(kk in 1:(nr * nc)) {
      mmk <- lapply(orders[[i]], derivatives)
      derivs[kk] <- eval(mmk[[kk]]) 
    }
      mt <- matrix(derivs, nrow = i)
      mtCol <- ncol(mt)
      MHval <- apply(mt, 2, prod) * coeffs[[i]]
      MH[i] <- sum(MHval)
    }
    return(MH)
  } 
  Lambda <- function(n, expr = expression(u^2), wrt1 = "a", 
    wrt2 = "u", aVal = 2, bVal = 1, xVal = 2, varU = TRUE, 
      lower = TRUE) {
    derivatives <- function(order) { 
      d <- Deriv(expr = expr, wrt = wrt2, order = order)  
      return(d)                               
    } 
    funCoeff <- function(n) {
      outCoeff <- vector('list', n)
      for (i in 1 : n) {
        k <- (n + i)
        outCoeff[[i]] <- Bell(n)[[n + i]]
      }
       out <- outCoeff 
       return(out)
    }
    funDerivOrder <- function(n) {
      outDerivOrder <- vector('list', n)
      nrowVec <- numeric(n)
      outCoeff <- vector('list', n)
      for (i in 1 : n) {
        outDerivOrder[[i]] <- Bell(n)[[i]]
        nrowVec[i] <- nrow(outDerivOrder[[i]]) 
      }
      out <- outDerivOrder
      return(out)
    }
    if(wrt1 == "a") {
      val <- aVal
      expr1 <-  expression(x * sqrt(a))
    } else {
      val <- bVal
      expr1 <- expression(sqrt(b)/x)
    }  
    if(lower == TRUE && varU == TRUE) {
      u <- sqrt(b)/x - x * sqrt(a)
      if(wrt1 == "a") {
        val <- aVal
        expr1 <-  expression(- x * sqrt(a))
      } else {
        val <- bVal
        expr1 <-  expression(sqrt(b)/x)
      }  
    } 
    if(lower == FALSE && varU == TRUE) {
      u <- x * sqrt(a) - sqrt(b)/x
      if(wrt1 == "a") {
      val <- aVal
      expr1 <-  expression(x * sqrt(a))
      } else {
      val <- bVal
      expr1 <- expression(- sqrt(b)/x)
      }  
    } 
    MH <- numeric(n)
    orders <- funDerivOrder(n)
    coeffs <- funCoeff(n) 
    for(i in 1:n) {
      nr <- nrow(orders[[i]])
      nc <- ncol(orders[[i]])
      derivs <- numeric(nr * nc)
      for(kk in 1:(nr * nc)) {
        mmk <- lapply(orders[[i]], derivatives)
        derivs[kk] <- eval(mmk[[kk]])          
      }
      mt <- matrix(derivs, nrow = i)
      mtCol <- ncol(mt)
      MHval <- apply(mt, 2, prod) * coeffs[[i]]
      MH[i] <- sum(MHval)
    }
    return(MH)
  }   
  Deriv <- function(expr, wrt, order = 1) {
    if(order == 0) return(expr)
    if(order == 1) D(expr, wrt)
    else Deriv(D(expr, wrt), wrt, order - 1)
  }   
  Erfc <- function(x) {
    2 * pnorm(-x * sqrt(2), lower.tail = TRUE)   
  }
  if(wrt1 == "a") {
    val <- aVal
    expr1 <-  expression(x * sqrt(a))
  } else {
    val <- bVal
    expr1 <- expression(sqrt(b)/x)
  }  
  if(lower == TRUE && varU == TRUE){
    u <- sqrt(b)/x - x * sqrt(a)
    if(wrt1 == "a") {
      val <- aVal
      expr1 <-  expression(- x * sqrt(a))
    } else {
      val <- bVal
      expr1 <-  expression(sqrt(b)/x)
    }  
  } 
  if(lower == FALSE && varU == TRUE){
    u <- x * sqrt(a) - sqrt(b)/x
    if(wrt1 == "a") {
    val <- aVal
    expr1 <-  expression(x * sqrt(a))
    } else {
    val <- bVal
    expr1 <- expression(- sqrt(b)/x)
    }  
  } 
  if(n >= 1) {    
    Delta <- Delta(n, expr = expr1, wrt = wrt1, val = val, 
      xVal = xVal)
    if (n >= 2){
      sum1 <- numeric(n - 1)
      for(beta in 2:n){
        LambdaOrder <- beta - 1
        sum2 <- numeric(LambdaOrder)
        for(i in 1:LambdaOrder){      
         sum2[i] <- (-1)^i * Lambda(LambdaOrder, expr = expr2, 
          wrt1 = wrt1, wrt2 = wrt2, aVal = aVal, bVal = bVal, 
          xVal = xVal, varU = varU, lower = lower)[i]
         sum1[beta] <- sum((sum2) * Delta[beta])
        }  
        re <- (- 2 * exp(-u^2))/sqrt(pi) * (Delta[1] + sum(sum1))
      }
    } else {  
      re <- (-2 * exp(-u^2))/sqrt(pi) * Delta[1]
    }
  }  
  if(n == 0) re <- Erfc(u)
  return(re)
} 
 
CalIncLapInt <- function(lambda, a = 1, b = 1, x = 1, lower = TRUE, 
  bit = 200) {
  stopifnot(a > 0, b >= 0, x > 0) 
  Deriv <- function(expr, wrt, order = 1) {
    if(order == 0) return(expr)
    if(order == 1) D(expr, wrt)
    else Deriv(D(expr, wrt), wrt, order - 1)
  }   
  derivRaExp <- function(k, expr1 = expression(exp(-2 * sqrt(a * b))), 
    wrt = "a", expr2 = expression(1/sqrt(a)), aVal = 1, bVal = 3) {   
    a <- aVal
    b <- bVal
    re <- numeric(k+1)  
    for(ss in 0:k) {
      d1 <- eval(Deriv(expr = expr1, wrt = wrt, order = k - ss))
      d2 <- eval(Deriv(expr = expr2, wrt = wrt, order = ss))
      re[ss + 1] <- choose(k, ss) * d1 * d2
    }  
    out <- sum(re)
    return(out)
  } 
  if(lambda < 0) {
    j <- -(2 * lambda + 1)/2
    wrt1 <- "a"
    varU1 <- TRUE     
    varU2 <- FALSE 
    expr <- expression(1/sqrt(a))
    if(lower == TRUE) {     
      ind1 <- 1         
    } else {            
      ind1 <- -1        
    }
  }
  if(lambda > 0) {
    j <- (2 * lambda - 1)/2 
    wrt1 <- "b"
    varU1 <- TRUE     
    varU2 <- FALSE 
    expr <- expression(1/sqrt(b))
    if(lower == TRUE) {     
      ind1 <- -1         
    } else {            
      ind1 <- 1        
    }
  }   
  if(j > 19)stop("Computational limit reached, abs of 
  lambda must be < 19/2")                  
  if(!is.whole(j))stop("lambda must be half an odd integer")
  ibk <- numeric(j + 1)
  ibk <- mpfr(ibk, bit)
    for(k in 0:j) {
      n1 <- mpfr(derivErfc(j - k, wrt1 = wrt1, 
            expr2 = expression(u^2), wrt2 = "u", aVal = a, 
            bVal = b, xVal = x, varU = varU1, lower = lower), bit)
      n3 <- mpfr(derivErfc(j - k, wrt1 = wrt1, 
            expr2 = expression(u^2), wrt2 = "u", aVal = a, 
            bVal = b, xVal = x, varU = varU2, lower = lower), bit)   
      if(b != 0) {
        n2 <- mpfr(derivRaExp(k, 
              expr1 = expression(exp(-2 * sqrt(a * b))), 
              expr2 = expr, wrt = wrt1, aVal = a, bVal = b), bit) 
        n4 <- mpfr(derivRaExp(k, 
              expr1 = expression(exp(2 * sqrt(a * b))), 
              expr2 = expr, wrt = wrt1, aVal = a, bVal = b), bit)  
        ibk[k + 1] <- chooseMpfr(j, k) * (n1 * n2 - ind1 * n3 * n4)  
      } 
      if(b == 0) {
        n5 <-  eval(Deriv(expr = expression(1/sqrt(a)), wrt = wrt1, 
          order = k))
        ibk[k + 1] <- chooseMpfr(j, k) * (n1 * n5 - ind1 * n3 * n5)  
    }  
  }  
  cons <- (-1)^j * sqrt(pi)/4 
  out <- as.numeric(mpfr(cons, bit) * sum(ibk))
  return(out)
}
  
pgig <- function(q, lambda, chi, psi, lower.tail = TRUE, bit = 200) {
  stopifnot(q > 0, chi > 0, psi > 0)
  z <- sqrt(chi * psi)    
  x <- sqrt(1/(2 * q * psi))
  if(lower.tail == FALSE) {
    CalIncLapInt(lambda = lambda, a = z^2, b = 1/4, x = x, lower = TRUE, 
      bit = bit)/((2 * z)^(lambda) * 
      besselK(z, nu = lambda, expon.scaled = FALSE))
  } else {
    CalIncLapInt(lambda = lambda, a = z^2, b = 1/4, x = x, lower = FALSE, 
      bit = bit)/((2 * z)^(lambda) * besselK(z, nu = lambda, 
      expon.scaled = FALSE))
  }                       
}

gamma_inc_err <- function(x, lambda, bit, lower = FALSE) {
  stopifnot(x > 0, lambda > 0)
  if(lower == TRUE){
    2 * CalIncLapInt(lambda = -lambda, a = 1, b = 0, x = sqrt(x), 
      lower = TRUE, bit = bit)
  } else {
    2 * CalIncLapInt(lambda = -lambda, a = 1, b = 0, x = sqrt(x), 
      lower = FALSE, bit = bit)
  }
}

besselK_inc_err <- function(x, z, lambda, bit, lower = FALSE) {
  if(lower == TRUE){
    CalIncLapInt(lambda = lambda, a = z^2, b = 1/4, x = x, 
      lower = TRUE, bit = bit)/(2 * z)^(lambda)
  } else {
    CalIncLapInt(lambda = lambda, a = z^2, b = 1/4, x = x, 
      lower = FALSE, bit = bit)/(2 * z)^(lambda)
  }
}

besselK_inc_ite <- function(x, y, lambda, traceIBF = TRUE, 
  epsilon = 0.95, nmax = 120) {
  Aki <- function(nmax, lambda) {
  A <- matrix(rep(0, (nmax + 1)^2), ncol = nmax + 1)
  A[1,1] <- 1
    for(k in 2:(nmax+1)) {
      for(i in 2:k) {
        A[k, i] <- (-lambda + i + k - 3) * A[k-1, i] + A[k-1, i-1]
      }
        A[k, k] <- 1
        A[k, 1] <- (-lambda + k - 2) * A[k-1,1]
      }
    return(A)
  }
  PT <- function(lambda) {
    Cnp <- numeric(0.5 * (lambda + 1) * (lambda + 2))
    Cnp[1:min(3, 0.5 *(lambda + 1) * (lambda + 2))] <- 1
    if (lambda > 1) {
      for(n in 2:lambda) {
        Cnp[n * (n + 1)/2 + 1] <- 1
        Cnp[n * (n + 1)/2 + n + 1] <- 1
        for(np in 1:(n - 1) ) {
          Cnp[n * (n+1)/2 + np + 1] <-
          Cnp[n * (n-1)/2 + np] + Cnp[n * (n-1)/2 + np + 1]
        }
      }
    }
  return(Cnp)
  }
  DENOM <- function(n, x, y ,lambda, An, nmax, Cnp){
    GN <- 0
      for (j in 0:n){
        terme <- 0
        for (i in 0:j){
          terme <- terme + An[j + 1, i + 1] * x^i
        }
        GN <- GN + Cnp[n * (n + 1)/2 + j + 1] * (-1/y)^j * terme
      }
      GN <- GN  * (-x * y)^n * x^(lambda + 1) * exp(x + y)
      return(GN)
  }
  NUM <- function(n, x, y, lambda, Am, An, nmax, Cnp, GM, GN){
    GM <- 0
    for (j in 1:n){
      terme <- 0
      for (k in 0:(j - 1)){
        termepr <- 0
        for (i in 0:k){
          termepr <- termepr + Am[k + 1, i + 1] * (-x)^i
        }
          terme <- terme + termepr * Cnp[j * (j - 1)/2 + k + 1] * (1/y)^k
        }
        GM <- GM + Cnp[n * (n + 1)/2 + j + 1] * (x * y)^j * GN[n - j + 1] * 
          terme
      }
      GM <- GM  * exp(-(x + y)) * x^(-lambda)/y
      return(GM)
    }
    nmax <- nmax
    tol <- .Machine$double.eps^epsilon
    if (epsilon >= 1) stop("epsilon should be less than 1")
    if (epsilon <= 0.5) stop("epsilon should be greater than 0.3")
    Am <- matrix(rep(0, (nmax + 1)^2), ncol = nmax + 1)
    An <- matrix(rep(0, (nmax + 1)^2), ncol = nmax + 1)
    Cnp <- numeric((nmax + 1) * (nmax + 2)/2)
    G <- numeric(nmax)
    GN <- numeric(nmax + 1)
    GM <- numeric(nmax)
    Cnp <- PT(nmax)  
    BK <- besselK(2 * sqrt(x * y), lambda)  
    Delta <- Omega <- Lambda<- 1
    fmax <- 1.797693e+308 
    fmin <- 2.225074e-308
    epsilon <- 2.220446e-16  
    if (BK < epsilon^0.90) Delta <- 0 
    if(x >= y) {                                                                 
      Am <- Aki(nmax, lambda - 1)                                             
      An <- Aki(nmax, -lambda -1)                                             
      GN[1] <- DENOM(0, x, y, lambda, An, nmax, Cnp)                         
      GN[2] <- DENOM(1, x, y, lambda, An, nmax, Cnp)                         
      GM[1] <- NUM(1, x, y, lambda, Am, An, nmax, Cnp, GM, GN)               
      G[1] <- x^lambda * GM[1]/GN[2]    
      for(n in (2:nmax)) {                                                  
        GN[n + 1] <- DENOM(n, x, y, lambda, An, nmax, Cnp)     
        GM[n] <- NUM(n, x, y, lambda, Am, An, nmax, Cnp, GM, GN)    
        G[n] <- x^lambda * GM[n]/GN[n + 1]                            
        nume <- x^lambda * GM[n]                                      
        deno <- GN[n + 1]
        if (is.infinite(nume)) {
            stop(paste("Infinity occurs in the numerator of G_n at n = ", 
              n, sep = ""))
        } 
        if (is.nan(nume)) {
            stop(paste("NaN occurs in the numerator of G_n at n = ", 
              n, sep = ""))
        }       
        if (is.infinite(deno)) {
            stop(paste("Infinity occurs in the denominator of G_n at n = ", 
              n, sep = ""))
        } 
        if (is.nan(deno)) {
            stop(paste("NaN occurs in the denominator of G_n at n = ", 
              n, sep = ""))
        }                                                        
        if (traceIBF == TRUE) {                                           
            cat("Iteration n    ", n, "\n")                               
            cat("G[n] numerator     :", nume, "\n")                  
            cat("G[n] denominator   :", deno, "\n")                  
            cat("Error              :", abs(G[n] - G[n-1]), "\n")    
        }
        if (abs(G[n] - G[n-1]) <= tol) {                                  
            break                                                         
        }                                                 
      }                                                                     
    }                                                                         
    if(y > x) {
      u <- sqrt(x * y)
      BK <- besselK(2 * u, lambda)  
      Am <- Aki(nmax, -lambda - 1)
      An <- Aki(nmax, lambda -1)
      GN[1] <- DENOM(0, y, x, -lambda, An, nmax, Cnp)
      GN[2] <- DENOM(1, y, x, -lambda, An, nmax, Cnp)
      GM[1] <- NUM(1, y, x, -lambda, Am, An, nmax, Cnp, GM, GN)
      G[1] <- y^(-lambda) * GM[2]/GN[1]
      for(n in (2:nmax)) {
        GN[n + 1] <- DENOM(n, y, x, -lambda, An, nmax, Cnp)
         GM[n] <- NUM(n, y, x, -lambda, Am, An, nmax, Cnp, GM, GN)
         G[n] <- y^(-lambda) * GM[n]/GN[n + 1]
         nume <- y^(-lambda) * GM[n]
         deno <- GN[n + 1]
         if (is.infinite(nume)) {
             stop(paste("Infinity occurs in the numerator of G_n at n = ", 
              n, sep = ""))
         } 
         if (is.nan(nume)) {
             stop(paste("NaN occurs in the numerator of G_n at n = ", 
              n, sep = ""))
         }       
         if (is.infinite(deno)) {
             stop(paste("Infinity occurs in the denominator of G_n at n = ", 
              n, sep = ""))
         } 
         if (is.nan(deno)) {
             stop(paste("NaN occurs in the denominator of G_n at n = ", 
              n, sep = ""))
         }                   
         if (traceIBF == TRUE) {                                                    
             cat("Iteration n    ", n, "\n")                                   
             cat("G[n] numerator     :", nume, "\n")                           
             cat("G[n] denominator   :", deno, "\n")
             cat("Error              :", abs(G[n] - G[n-1]), "\n") 
         }
         if (abs(G[n] - G[n-1]) <= tol) {   
          G[n] <- 2 * ((x/y)^(lambda/2)) * BK - G[n]                               
          break                                                         
        }            
      }  
    }
    if(Delta == 0) warning("Loss of accuracy because the algorithm 
      stops prematurely")
    inBK <- G[n]
    return(inBK)
}    

besselK_inc_clo <- function (x, z, lambda, lower = FALSE, 
  expon.scaled = FALSE) {
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)
    abs(x - round(x)) < tol
  if (z >= 705 & expon.scaled == FALSE)
    warning("Underflow occurs, expon.scaled should be TRUE")
  if (z == 0)
    stop ("z should be positive")
  if ( x < 0)
    stop("x should be greater than zero")
  if (lambda < 0)
    stop("lambda should be positive")
  if (!is.wholenumber(lambda + 1/2))
    stop("lambda should be half of an odd positive integer")
  s <- lambda + 1/2
  Az <- sqrt(pi/(2 * z)) * exp(-z)
  if (expon.scaled == TRUE)
    Az <- sqrt(pi/(2 * z))
  r <- seq(0, s-1, length.out = s)
  evalbk <- rep(0, length(r))
  for(i in 1:length(evalbk)) {
    alpha <- lambda + r[i] + 1/2
    beta <- factorial(lambda - r[i]- 1/2)
    varphi <- seq(0, alpha - 1, length.out = alpha)
    binoCoef <- choose(alpha-1, varphi)
    binoVal <- numeric(length(binoCoef))
    for (j in 1:length(binoVal)) {
      binoVal[j] <- binoCoef[j] * x^(-varphi[j]) * 
      factorial(varphi[j]) 
    }
    if (lower== FALSE) {
      evalbk[i] <- 1/(beta*factorial(r[i]))*
        (1/(2 * z))^(r[i]) * x^(alpha - 1) * exp(-x) * sum(binoVal)
    } else {
      evalbk[i] <- 1/(beta*factorial(r[i]))*
        (1/(2 * z))^(r[i]) * (factorial(alpha - 1) - 
        x^(alpha - 1)*exp(-x) * sum(binoVal))
    }   
  }
  Az * sum(evalbk)
}
