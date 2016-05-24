perm.f.test <-
function(response, treatment=NULL, num.sim=20000)  {
  # A permutation F-test is performed, and a one-way analysis of variance F-test is performed.
  # Note: Only approximated p-values are calculated for the permutation F-test.
  # 'response': Numeric vector of responses if treatment is not \code{NULL}.
  #             If \code{treatment} is \code{NULL}, then \code{response} must be an N by 2 matrix,
  #             such that the first column represents response and the second column represents treatment.
  # 'treatment': Vector of treatments, which need not be numerical.
  #              code{treatment} should be set to \code{NULL} if \code{response} is an N by 2 matrix.
  # 'num.sim': The number of simulations performed.
  #            If \code{num.sim} is smaller than one, then the permutation test is not performed.
  # Example:  perm.f.test( c( 3,6,5,2,14,7,9,15,11,13,12 ), rep( 1:3, c(4,4,3) ) )
  if (is.null(treatment))   {treatment <- response[,2] ;  response <- as.numeric(response[,1])}
  if (length(response)!=length(treatment))
    stop("'response' and 'treatment' must have the same length.")
  if ( !is.numeric(num.sim) )  stop("'num.sim' must be numeric.")
  treatment <- factor(treatment) ;  list1 <- c("One-way ANOVA", summary(aov(response~treatment)))
  if (num.sim>=1) {
    all.treatments <- union(treatment, NULL); k <- length(all.treatments); ni <- NULL
    for (i in 1:k) { ni <- c(ni, sum(treatment==all.treatments[i])) }
    one.way.anova <- function(treatment, response, all.treatments, ni, k) {
         # Returns $\sum_{i=1}^k n_i{\bar Y}_i^2$, which is monotone increasing with          
         #   the F-test statistic associated with one-way ANOVA; see Higgins, 2004, p. 84.
         yi.mean <- NULL
         for (i in 1:k) {yi.mean <- c(yi.mean, mean(response[treatment==all.treatments[i]])) }
         sum(ni*yi.mean^2)   }
    test.stat0 <- one.way.anova(treatment, response, all.treatments, ni, k) ;   test.stat <- NULL
    for (i in 1:num.sim) {
      test.stat <- c(test.stat, one.way.anova(treatment, sample(response), all.treatments, ni, k)) }
    list2 <- c( "The p-value from the permutation F-test based on",
                paste( num.sim, "simulations is", mean(test.stat>=test.stat0) ) )
    list1 <- c(list1, list2) }
  return(list1) 
}
