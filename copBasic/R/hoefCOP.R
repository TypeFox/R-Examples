"hoefCOP" <-
function(cop=NULL, para=NULL, p=2, as.sample=FALSE, sample.as.prob=TRUE,
                              brute=FALSE, delta=0.002, ...) {

    if(as.sample) {
      if(is.null(para)) {
         warning("Sample Hoeffding's Phi desired but para is NULL, ",
                 "returning NULL")
         return(NULL)
      }
      if(length(names(para)) != 2) {
        warning("para argument must be data.frame having only two columns, ",
                "returning NULL")
        return(NULL)
      }
      if(as.sample == -1) message("Sample Hoeffding's Phi---CPU intensive!")

      # Gaisser et al. (2010, eq. 3): A multivariate extension version of
      #  Hoeffding's Phi-Square: J. Multivariate Analysis: v. 101, pp. 2571--2586
      "invHd" <- function(d) {
         if(d <= 0) stop("dimension can not be 0 or negative!")
         HdA <- 2/((d+1)*(d+2)) + (1/3)^d
         HdB <- (1/2^d) * factorial(d) /
                   cumprod(sapply(0:d, function(i) { i + 1/2 }))[d+1] # d + 1 !!!!
         return(HdA - HdB)
         # The value 1/invHD(2) must be 90.
      }

      # The equation by those authors is based on probabilities for U. There
      # is a chance that users might feed non probability values so in that
      # case, rank and convert back to a simple probability.
      U <- as.matrix(para); n <- length(U[,1]); d <- length(U[1,])
      if(! sample.as.prob) for(i in 1:d) U[,i] <- rank(U[,i])/n
      A <- sum(sapply(1:n, function(j) {
             sum(sapply(1:n, function(k) {
             cumprod(sapply(1:d, function(i) { 1 - max(c(U[j,i], U[k,i])) }))[d]
             } ))
           } ))

      B <- sum(sapply(1:n, function(j) {
                           cumprod(sapply(1:d, function(i) { 1 - U[j,i]^2 }))[d]
           } ))

      samPHIsq <- ((1/n)^2) * A - ((2/n)*(1/2)^d)*B + (1/3)^d
      hd <- 1/invHd(d)
      # Note that Hoeffding's Phi (not Phi-Square is being returned)
      return(sqrt(hd*samPHIsq))
   }

   if(is.null(cop)) {
      warning("must have copula argument specified, returning NULL")
      return(NULL)
   }

   if(p == 1) {
      kp <- 12
   } else if(p == 2) {
      kp <- 90
   } else if(p == 3) {
      kp <- 560
   } else if(p == 4) {
      kp <- 3150
   } else if(p == 5) {
      kp <- 16600
   } else {
      kp <- exp(lgamma(2*p+3) - 2*lgamma(p + 1))/2 # Nelsen (2006, p. 213)
   }
   if(brute) {
      us <- vs <- seq(.Machine$double.eps, 1-.Machine$double.eps, delta)
      hoef <- sum(sapply(us, function(u) {
                 sum(sapply(vs, function(v) {
                    abs(COP(u,v, cop=cop, para=para, ...) - u*v)^p }))
               }))
      hoef <- (hoef*delta^2*kp)^(1/p)
      return(hoef)
   }

   myint <- NULL
   try(myint <- integrate(function(u) {
            sapply(u,function(u) {
                      integrate(function(v) {
                          abs(COP( u, v, cop=cop, para=para, ...) - u*v)^p},
                      0, 1)$value
            })}, 0, 1) )
   if(is.null(myint)) {
      return(NA)
   } else if(myint$value == 0 & myint$abs.error == 0) {
      warning("integrate() returned zero with zero absolute error, ",
              "p is likely too large or the copula is Independence")
      return(NA)
   } else {
      return((myint$value*kp)^(1/p))
   }
}

"LpCOP" <-
function(cop=NULL, para=NULL, p=2, brute=FALSE, delta=0.002, ...) {
    return(hoefCOP(cop=cop, para=para, p=p, brute=brute, delta=delta, ...))
}


"LpCOPradsym" <-
function(cop=NULL, para=NULL, p=2, brute=FALSE, delta=0.002, ...) {
   if(brute) {
      sum <- 0
      us <- vs <- seq(.Machine$double.eps, 1-.Machine$double.eps, delta)
      hoef <- sum(sapply(us, function(u) {
                 sum(sapply(vs, function(v) {
                    abs(COP(u,v, cop=cop, para=para, ...) -
                    surCOP(u,v, cop=cop, para=para, ...))^p }))
              }))
      kp <- 1
      Lp <- (sum*hoef^2*kp)^(1/p)
      return(Lp)
   }

   myint <- NULL
   try(myint <- integrate(function(u) {
            sapply(u,function(u) {
                      integrate(function(v) {
                          abs(COP( u, v, cop=cop, para=para, ...) -
                           surCOP(u,v, cop=cop, para=para, ...))^p},
                      0, 1)$value
            })}, 0, 1) )
   if(is.null(myint)) {
      return(NA)
   } else {
      kp <- 1
      return((myint$value*kp)^(1/p))
   }
}


"LpCOPpermsym" <-
function(cop=NULL, para=NULL, p=2, brute=FALSE, delta=0.002, ...) {
   if(brute) {
      us <- vs <- seq(.Machine$double.eps, 1-.Machine$double.eps, delta)
      hoef <- sum(sapply(us, function(u) {
                 sum(sapply(vs, function(v) {
                    abs(COP(u,v, cop=cop, para=para, ...) -
                        COP(v,u, cop=cop, para=para, ...))^p }))
              }))
      kp <- 1
      Lp <- (hoef*delta^2*kp)^(1/p)
      return(Lp)
   }

   myint <- NULL
   try(myint <- integrate(function(u) {
            sapply(u,function(u) {
                      integrate(function(v) {
                          abs(COP(u,v, cop=cop, para=para, ...) -
                              COP(v,u, cop=cop, para=para, ...))^p},
                      0, 1)$value
            })}, 0, 1) )
   if(is.null(myint)) {
      return(NA)
   } else {
      kp <- 1
      return((myint$value*kp)^(1/p))
   }
}

