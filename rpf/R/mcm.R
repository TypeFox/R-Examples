##' Create a multiple-choice response model
##'
##' WARNING: This model is mostly not implemented.
##' 
##' This function instantiates a multiple-choice response
##' model.
##' 
##' @param outcomes the number of possible outcomes
##' @param numChoices the number of choices available
##' @param factors the number of factors
##' @return an item model
##' @export
##' @author Jonathan Weeks <weeksjp@@gmail.com>
rpf.mcm <- function(outcomes=2, numChoices=5, factors=1) {
  stop("Not implemented")
  guess.weight <- 20
  guessing <- (1/numChoices)
  new("rpf.mdim.mcm",
      outcomes=outcomes, factors=factors,
      a.prior.sdlog=.5,
      c.prior.alpha=guess.weight*guessing+1,
      c.prior.beta=guess.weight*(1-guessing)+1)
}

### mdim

setMethod("rpf.prob", signature(m="rpf.mdim.mcm", param="numeric",
                                theta="matrix"),
          function(m, param, theta) {
            den <- NULL
            
            a1 <- param[1:(m@factors*m@outcomes)]
            b1 <- param[(m@factors*m@outcomes+1):
                        ((m@factors+1)*m@outcomes)]
            c1 <- param[-1:-((m@factors+1)*m@outcomes)]
            
            ##   Compute the denominator
            for (k in 1:m@outcomes) {
              tmp <- (k-1)*m@factors
              tmp1 <- tmp+m@factors
              d <- exp((theta %*% a1[(tmp+1):tmp1])+b1[k])
              den <- cbind(den, d)
            }
            den <- apply(den,1,sum)

            numPersons <- dim(theta)[1]

            p <- array(dim=c(numPersons, m@outcomes))

            for (k in 1:m@outcomes) {
              tmp <- (k-1)*m@factors
              tmp1 <- tmp+m@factors
              if (k==1) {
                cp <- (exp((theta %*% a1[(tmp+1):tmp1])+b1[k]))/den
              } else {
                cp <- (exp((theta %*% a1[(tmp+1):tmp1])+b1[k])+c1[k-1]*(exp((theta %*% a1[1:m@factors])+b1[1])))/den
              }
              p[,k] <- cp
            }
            return(p)
          })

setMethod("rpf.rparam", signature(m="rpf.mdim.mcm"),
          function(m) {
              a <- rlnorm(m@outcomes * m@factors,
                          meanlog=0, sdlog=m@a.prior.sdlog)
              b <- sort(rnorm(m@outcomes))
              c <- rbeta(m@outcomes-1, shape1=m@c.prior.alpha-2,
                         shape2=m@c.prior.beta-2)
              c(a=a,b=b,c=c)
          })
