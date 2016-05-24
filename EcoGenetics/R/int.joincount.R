#' Join-count statistic, internal.
#' 
#' @param Z Vector, matrix or data frame.
#' @param con Connection network.
#' @param ncod number of elements coding each category (e.g., if x ncod =1, if xx, 
#' ncod = 2, and so on).
#' @param nsim Number of Monte-Carlo simulations. 
#' @param test If test = "bootstrap", the program generates a bootstrap 
#' resampling and the associated confidence intervals of the null hypothesis.
#'  If test = "permutation" (default) a permutation test is made and the p value 
#'  is calculated. 	 
#' @param alternative The alternative hypothesis. If "auto" is selected (default) the
#' program determines the hypothesis by difference between the median of the simulations
#' and the observed value. Other options are: "two.sided", "greater" and "less".
#' if test == cross, for the first interval (d== 0) the p and CI are computed with cor.test.
#' @param adjust.n Should be adjusted the number of individuals? (warning, this would
#' change variances)
#' @param adjust Method for multiple correction of P-values 
#' passed to \code{\link[stats]{p.adjust}}.
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' 
#' @keywords internal

int.joincount <- function(Z, con, ncod, ploidy, nsim,
                          alternative, test = "permutation", 
                          adjust.n = FALSE, adjust) {
   
  con <- int.check.con(con)
  con <- as.vector(con)
  
  if(test == "permutation") {
    replace <- FALSE
  } else {
    replace <- TRUE
  }
  
  
  #internal function 
  
  jcfun <- function(input) {
    
    outmat <- outer(input, input, FUN = "paste", sep = "")
    outmat <- as.matrix((aue.sort(outmat, ploidy = 2)))  #symmetric matrix
    outmat <- as.factor(outmat)

    temp <- list()
    for(i in seq(along = levels(outmat))) {
      temp[[i]] <- as.integer(outmat == levels(outmat)[i])
    }
    
    sapply(temp, function(x) sum(x * con) / 2) 
  }
  
  obs <- jcfun(Z)
  
  #simulated datasets. permuting rows and columns, mantaining structure
  monte.c <- matrix(0, nrow = nsim, ncol = length(obs))
  for(i in 1:nsim) {
    samp <- sample(length(Z), replace = replace)
    outsamp <- Z[samp]
    monte.c[i, ] <- jcfun(outsamp)
  }
  
  ran <- int.random.test(repsim = monte.c, obs = obs, 
                         nsim = nsim, test = test,
                         alternative = alternative,
                         adjust = adjust)
  
  
  #labeling rows
  outmat <- outer(Z, Z, FUN = "paste", sep = "")
  rownames(ran)<- levels(as.factor(aue.sort(outmat, ploidy = 2))) 
  
  res <- list("analysis" = "Join-count", 
              "nsim" = nsim,
               "results" = ran)
  
  res
  
}
