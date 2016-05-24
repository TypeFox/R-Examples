#$Author: sinnwell $
#$Date: 2013/01/14 19:32:41 $
#$Header: /projects/genetics/cvs/cvsroot/haplo.stats/R/haplo.binomial.R,v 1.1 2013/01/14 19:32:41 sinnwell Exp $
#$Locker:  $
#$Log: haplo.binomial.R,v $
#Revision 1.1  2013/01/14 19:32:41  sinnwell
#add haplo.binomial, no need glm.fit.nowarn
#


## Purpose: Extracted from R 2.15.0 to "fix" the init expression to not give warnings
##          for binomial fits beause of our non-integer weights, which are
##          posterior prob of haplotypes for a subject


haplo.binomial <- function (link = "logit") 
{
    save <- binomial()
    save$initialize <- expression({
        if (NCOL(y) == 1) {
            if (is.factor(y)) y <- y != levels(y)[1L]
            n <- rep.int(1, nobs)
            y[weights == 0] <- 0
            ## The only line that is different from binomial() in R 2.15.0 is
            ## disabling the lines that warn of non-integer counts
            ## if (any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1")
            mustart <- (weights * y + 0.5)/(weights + 1)
            m <- weights * y
            ##if (any(abs(m - round(m)) > 0.001))
              ## warning("non-integer #successes in a binomial glm!")
        } else if (NCOL(y) == 2) {
            ## if (any(abs(y - round(y)) > 0.001))
              ## warning("non-integer counts in a binomial glm!")
            n <- y[, 1] + y[, 2]
            y <- ifelse(n == 0, 0, y[, 1]/n)
            weights <- weights * n
            mustart <- (n * y + 0.5)/(n + 1)
        } else stop("for the binomial family, y must be a vector of 0 and 1's\n", 
"or a 2 column matrix where col 1 is no. successes and col 2 is no. failures")
    })

    return(save)
    
}
