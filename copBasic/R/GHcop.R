# Gumbel-Hougaard Copula
"GHcop" <-
function(u, v, para=NULL, tau=NULL, tau.big=0.985, cor=NULL, ...) {

  if(! is.null(cor)) tau <- cor

  # tau.big <- 0.985 # tests show at least a failure rate < 1/10000 in simCOP
  # if tau is kept this low, which is a really large tau.
  para.big <- 1/(1-tau.big)

  if(is.null(para)) {
     if(is.null(tau)) {
        if(    length(u)     ==  length(v) &
              (length(u) > 1  &  length(v) > 1)) {
           tau <- cor(u,v, method="kendall")
        } else {
           warning("Argument para (Alpha) is NULL and u an v have unequal length",
                   "that prevents Alpha estimation by Kendall's Tau, returning NULL")
           return(NULL)
        }
        para <- 1/(1-tau)
        if(para < 1) {
           warning("negative correlation can not be fit with this copula, ",
                   "returning NULL")
           return(NULL)
        }
        names(para) <- "theta"; names(tau) <- "tau"
        return(list(para=para, tau=tau))
     } else {
        if(tau < 0) {
           warning("Kendall's Tau is negative, not compatible with this copula, ",
                   "returning NULL")
           return(NULL)
        }
        para <- 1/(1-tau)
        names(para) <- "theta"; names(tau) <- "tau"
        return(list(para=para, tau=tau))
     }
  }
  if(length(para) == 1) {
     if(para[1] < 1 | para[1] > para.big ) {
        warning("Parameter Theta (para) is not in [1, ",para.big,
                "] (numerically too big), returning M(u,v)")
        return(M(u,v))
     }
  } else if(para[1] < 1) {
     warning("Parameter Theta (para[1]) is < 1, returning NULL")
     # Because the two pis control as well, identification of a numerical
     # upper bounds of Theta for the asymmetric GH is not possible
     return(NULL)
  }
  if(length(u) > 1 & length(v) > 1 & length(u) != length(v)) {
    warning("length u = ",length(u), " and length v = ",length(v))
    warning("longer object length is not a multiple of shorter object length, ",
            "no recycling")
    return(NA)
  }
  # The extra hassle of vectorization made here is to handle situations
  # in which nested integrals are used where uneven vectors can be passed
  if(length(u) == 1) {
     u <- rep(u, length(v))
  } else if(length(v) == 1) {
     v <- rep(v, length(u))
  }

  if(length(para) == 1) {
     cop <- ((-log(u))^para + (-log(v))^para)^(1/para)
     return(exp(-cop))
  } else if(length(para) == 3) {
     pi2 <- para[2]; pi3 <- para[3]; di <- 1/para[1]
     x <- -log(u); y <- -log(v)
     cop <- exp(-(((pi2*x)^para[1] + (pi3*y)^para[1])^di + (1-pi2)*x + (1-pi3)*y))
     cop[is.nan(cop)] <- 0
     return(cop)
  } else {
     stop("Should not be here in logic")
  }
}


