size.selection.bechhofer <- function(beta, a, delta, sigma=NA)
{
  S <- as.matrix(pdIdent(diag(a-1)))
  S <- 0.5*(S+matrix(1,a-1,a-1))
  if (!is.na(sigma))
    {q <- qmvnorm(1-beta,tail="lower.tail",corr=S)$quantile
     b <- sigma^2/(delta^2)
     n <- ceiling(2*b*q^2)
     return(n)
    }
  else
    {b <- 1/(delta^2)
     g <- function(n){
            df <- a*(n-1)
            df <- as.integer(df)
            n-2*b*(qmvt(1-beta,tail="lower.tail",
                   df=df,corr=S)$quantile)^2
     }
     n <- uniroot(g,c(1,10^6))$root
     return(ceiling(n))
    }
}


                                        # requires library(mvtnorm)
                                        # requires library(nlme)
