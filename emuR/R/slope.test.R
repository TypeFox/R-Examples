##' Slope Test
##' 
##' Tests whether the difference between two or more regression lines is
##' significant
##' 
##' 
##' @param ...  this function takes any number of two column matrices. The
##' first column is the y-data (in the case of locus equations, this is the
##' vowel onset) and the second column is the x-data (in the case of locus
##' equations, vowel target).
##' @return The return value consists of the following componenets:
##' 
##' \item{separate}{ slope, intercept, r-squared, F-ratio, "d(egrees of)
##' f(reedom)" and "prob(ability that) line fits data" for the separate data
##' matrices entered. } \item{combined}{ F-ratio, "d(egrees of) f(reedom)", and
##' "Probability of them being DIFFERENT" for the slope and for the intercept
##' of the combined data.  } \item{x}{ the combined x-data for all the
##' matrices.  } \item{y}{ the combined y-data for all the matrices.  }
##' \item{mat}{ the category vectors for the combined data (consists of 1, 0
##' and -1).  } \item{numrows}{ the number of rows in each matrix.  }
##' \item{numcats}{ the sum number of matrices entered.
##' 
##' }
##' @seealso lm(), summary.lm(), pf()
##' @references see E. Pedhazur, Multiple Regression in Behavioral Research
##' p.436-450, 496-507.
##' @keywords misc
##' @export Slope.test
"Slope.test" <- function(...)
{
  ## compiled by Jonathan Harrington and Marija Tabain (October 1997)
  ## this function tests whether the intercepts and slopes
  ## of two or more (straight-line) regressions are significantly
  ## different
  ## matrices are to be compared on slope and intercept;
  ## arrange y in col 1, x in col 2 in each case.
  ## see E. Pedhazur, Multiple Regression in Behavioral Research
  ## p.436-450, 496-507. 
  Slope.sub <- function(...)
  {
    ## combine the matrices, and find out how many rows there are altogether
    omat <- NULL
    omat$numcats <- length(list(...))
    for(j in list(...)) {
      numrows <- nrow(j)
      omat$y <- c(omat$y, j[, 1])
      omat$x <- c(omat$x, j[, 2])
      omat$numrows <- c(omat$numrows, numrows)
    }
    ## set up category vectors of 1, 0, 0, .... -1
    vec <- rep(0, omat$numcats - 1)
    omat$mat <- NULL
    for(j in 1:length(vec)) {
      zeros <- vec
      zeros[j] <- 1
      zeros <- c(zeros, -1)
      zeros <- rep(zeros, omat$numrows)
      omat$mat <- cbind(omat$mat, zeros)
    }
    omat
  }
  
  ## main function begins here
  omat <- Slope.sub(...)	
  ## number of category vectors and the (1) continuous vector for intercept
  k1 <- omat$numcats	# the (1) continuous vector for intercept
  k2 <- 1	        # number of category vectors, product 
  # vectors and (1)continuous vector 
  
  ## for slope
  k3 <- 1 + ((omat$numcats - 1) * 2)	
  ## number of category vectors and (1) continuous vector for slope
  k4 <- omat$numcats	## length of y and of x
  N <- sum(omat$numrows)
  
  for(j in list(...)) {
    ## find the F-ratio, degrees of freedom, r-squared values, slope and intercept
    ## for the separate matrices
    firstvals <- summary.lm(lm(j[, 1] ~ j[, 2]))
    first.pf <- pf(firstvals$fstatistic[1], firstvals$fstatistic[2],
                   firstvals$fstatistic[3])
    first.out <- c(firstvals$r.squared, firstvals$fstatistic, 
                   first.pf, firstvals$coefficients[, 1])
    omat$separate <- rbind(omat$separate, first.out)
  }
  
  dimnames(omat$separate)[[2]] <- c("r-sq", "F ratio", "df", "df", 
                                    "prob. line fits data", "intercept", "slope")	
  
  ## multiply the category vectors by the x-values 
  prodvals <- omat$x * omat$mat
  z123 <- lm(omat$y ~ omat$x + omat$mat + prodvals)
  z12 <- lm(omat$y ~ omat$x + omat$mat)
  z2 <- lm(omat$y ~ omat$x)	## r-squared vals
  r12 <- summary.lm(z12)$r.squared
  r2 <- summary.lm(z2)$r.squared	## F-ratios
  fval.in.num <- (r12 - r2)/(k1 - k2)
  fval.in.den <- (1 - r12)/(N - k1 - 1)
  fratio.in <- fval.in.num/fval.in.den
  s123 <- summary.aov(z123)
  fratio.slope <- s123$"F Value"[3]	
  
  ## calculate probabilities and degrees of freedom
  prob.in <- pf(fratio.in, k1 - k2, N - k1 - 1)
  prob.slope <- pf(fratio.slope, k3 - k4, N - k3 - 1)
  first <- c(fratio.in, prob.in, k1 - k2, N - k1 - 1)
  second <- c(fratio.slope, prob.slope, k3 - k4, N - k3 - 1)
  outtemp <- rbind(first, second)
  col.lab <- c("intercept", "slope")
  row.lab <- c("F ratio", "Probability of them being DIFFERENT", "df", 
               "df")
  dimnames(outtemp) <- list(col.lab, row.lab)
  omat$combined <- outtemp
  omat
}
