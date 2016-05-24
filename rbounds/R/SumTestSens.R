SumTestSens <- function(T, q, n, m, Gamma) {
  
  ## Created by Devin Caughey on 26 February 2010
  ## Last modified on 4 April 2012

  ## SUMMARY
  ## The function 'SumTestSens' is an implementation of the method of
  ## sensitivity analysis for comparing two unmatched groups
  ## described in Section 4.6 of Paul Rosenbaum "Observational
  ## Studies" (2nd Ed., 2002).  It is designed to be used for sum
  ## statistics, such as Wilcoxon's rank sum  statistic.  An example
  ## of how this function may be used, taken from Section 4.6 of
  ## Rosenbaum (2002), is provided at the end of this code.

  ## ARGUMENTS
  ## 'T': observed value of the test statistic (e.g., the sum of
  ## the ranks of the responses of the treated group)
  ## 'q': vector of functions of the responses (e.g., their ranks;
  ## note that a higher rank corresponds to a higher response), sorted
  ## in  decreasing order (don't forget to do this).
  ## 'n': total number of observations
  ## 'm': number of treated observations
  ## 'Gamma': upper limit on the ratio of the a priori odds of
  ## treatment  assignment between the treated and control groups.

  ## RETURNS
  ## This function prints the upper bound of the normal approximation
  ## one-sided p-value for the test at the given value of Gamma. It
  ## also invisibly returns a list of intermediate statistics.
  
  K <- 0:n
  u <- matrix(nrow = n + 1, ncol = n)
  G <- Gamma
  g <- log(G)
  rho.i <- matrix(data = NA, nrow = n + 1, ncol = n)
  rho.ij <- array(data = NA, dim = c(n + 1, n, n))
  mu.T <- rep(NA, n + 1)
  var.T <- rep(NA, n + 1)
  sd.T <- rep(NA, n + 1)
  deviate <- rep(NA, n + 1) 

  ## Define function 'Z'.
    Z <- function(n, m, k, G) {
    Z.out <- 0
    if (m >= 0 & k >= 0) { 
      aa <- max(0, m + k - n):min(m, k)
      maa <- m - aa
      Z.out <- sum(choose(k, aa)*choose((n - k), maa)*(G^aa))
    }
    return(Z.out)
  }
  
  for(k in K) {
    ## Assign u[k + 1, ] k 0's followed by n - k 1's.
    u[k + 1, ] <- c(rep(1, k), rep(0, n - k))
    
    z_nmkG <- Z(n, m, k, G)
    for (i in 1:n) {
    ## Calculate unit i's probability of treatment
      rho.i[k + 1, i] <- exp(g * u[k + 1, i]) *
        Z(n - 1, m - 1, k - u[k + 1, i], G) / z_nmkG
      for (j in 1:n) {
        ## Calculate i and j's joint probability of treatment. 
        if (i == j) {
          rho.ij[k + 1, i, j] <- rho.i[k + 1, i]
        } else {
          rho.ij[k + 1, i, j] <-
            Z(n - 2, m - 2, k - u[k + 1, i] - u[k + 1, j], G) *
              exp(g*(u[k + 1, i] + u[k + 1, j])) / z_nmkG
        }
      } 
    }
    ## mean of T under the null
    mu.T[k + 1] <- q %*% rho.i[k + 1, ]
    ## standard deviation of T under the null
    var.T[k + 1] <- 0
    var.T.vec <- rep(NA, n)
    for(i in 1:n) {
      var.T.vec[i] <- sum(q[i]*q*(rho.ij[k + 1, i, ] -
                                  rho.i[k + 1, i]*rho.i[k + 1, ]))
    }
    var.T[k + 1] <- sum(var.T.vec)

    sd.T[k + 1] <- sqrt(var.T[k + 1])
    ## deviate
    deviate[k + 1] <- (T - mu.T[k + 1]) / sd.T[k + 1]
  }
  ## Main result
  minDeviate <- min(deviate)
  pValueUB <- pnorm(q = minDeviate,
                    mean = 0,
                    sd = 1,
                    lower.tail = FALSE)
  pValueUB.print <- ifelse(pValueUB < 0.0001,
                           "< 1e-04",
                           as.character(round(pValueUB, 4)))
  ## Collect output
  kMin <- K[which(deviate == min(deviate))]
  output <- list(pValueUB,
                 minDeviate,
                 deviate,
                 kMin,
                 T,
                 mu.T,
                 var.T,
                 sd.T,
                 rho.i,
                 rho.ij)
  names(output) <- c("pValueUB", "minDeviate", "deviate", "kMin", "T",
    "mu.T", "var.T", "sd.T", "rho.i", "rho.ij")
  ## Print p-value
  cat("For Gamma = ", G, 
      ", the upper-bound on the p-value of the sum test is: ", 
      pValueUB.print, ".\n", sep = "") 
  ## Return output invisibly
  invisible(output)
}

## Change FALSE to TRUE to test the following example from
## Rosenbaum (2002, p. 146)
if (FALSE) {
  mercury <- data.frame(matrix(c(1, 0, 2.7,    5.3,
                                 2, 0, 0.5,   15.0,
                                 3, 0, 0.0,   11.0,
                                 4, 0, 0.0,    5.8,
                                 5, 0, 5.0,   17.0,
                                 6, 0, 0.0,    7.0,
                                 7, 0, 0.0,    8.5,
                                 8, 0, 1.3,    9.4,
                                 9, 0, 0.0,    7.8,
                                10, 0, 1.8,   12.0,
                                11, 0, 0.0,    8.7,
                                12, 0, 0.0,    4.0,
                                13, 0, 1.0,    3.0,
                                14, 0, 1.8,   12.2,
                                15, 0, 0.0,    6.1,
                                16, 0, 3.1,   10.2,
                                17, 1, 0.7,  100.0,
                                18, 1, 4.6,   70.0,
                                19, 1, 0.0,  196.0,
                                20, 1, 1.7,   69.0,
                                21, 1, 5.2,  370.0,
                                22, 1, 0.0,  270.0,
                                23, 1, 5.0,  150.0,
                                24, 1, 9.5,   60.0,
                                25, 1, 2.0,  330.0,
                                26, 1, 3.0, 1100.0,
                                27, 1, 1.0,   40.0,
                                28, 1, 3.5,  100.0,
                                29, 1, 2.0,   70.0,
                                30, 1, 5.0,  150.0,
                                31, 1, 5.5,  200.0,
                                32, 1, 2.0,  304.0,
                                33, 1, 3.0,  236.0,
                                34, 1, 4.0,  178.0,
                                35, 1, 0.0,   41.0,
                                36, 1, 2.0,  120.0,
                                37, 1, 2.2,  330.0,
                                38, 1, 0.0,   62.0,
                                39, 1, 2.0,   12.8),
                               nrow = 39, ncol = 4, byrow = TRUE))
  colnames(mercury) <- c("ID", "Tr", "Pct.cu.cells", "Hg.in.blood")
  
  (T_test <- rank(mercury$Hg.in.blood) %*% mercury$Tr)
  (q_test <- sort(rank(mercury$Hg.in.blood), decreasing = TRUE))
  (n_test <- nrow(mercury))
  (m_test <- sum(mercury$Tr))

  ## Note: since this function uses exact rather than approximate
  ## formulas for the mean and variance of T, the p-values it
  ## calculates do not precisely match those in Rosenbaum (2002).
  testOut <- SumTestSens(T = T_test,
                         q = q_test,
                         n = n_test,
                         m = m_test,
                         Gamma = 1)
  ## Rosenbaum's value: < 1e-04

  testOut5 <- SumTestSens(T = T_test,
                          q = q_test,
                          n = n_test,
                          m = m_test,
                          Gamma = 5)
  ## > Rosenbaum's value: 0.0003

  testOut20 <- SumTestSens(T = T_test,
                           q = q_test,
                           n = n_test,
                           m = m_test,
                           Gamma = 20)
  ## > Rosenbaum's value: 0.0075
  
  testOut35 <- SumTestSens(T = T_test,
                           q = q_test,
                           n = n_test,
                           m = m_test,
                           Gamma = 35)
  ## > Rosenbaum's value: 0.0179

  ## Apply to vector
  sapply(c(1, 5, 20, 35), SumTestSens,
         T = T_test, q = q_test, n = n_test, m = m_test)
  
} ## end test


