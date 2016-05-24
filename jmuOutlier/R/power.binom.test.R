power.binom.test <-
function(n, alpha=0.05, alternative=c("two.sided", "less", "greater"), 
                                null.median, alt.pdist, ...) { 
  # Computes the power of the binomial test (i.e., sign test). 
  # This power function is based on binomial probabilities,
  #    not the normal approximation to the binomial distribution,
  #    and produces exact solutions.
  # 'n': The sample size.
  # 'alpha': The size of the test; i.e., P(type I error).
  # 'alternative': A character string specifying the alternative hypothesis, and
  #   must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}.  
  #   Only the initial letter needs to be specified.
  # 'null.median': The population median under the null hypothesis.
  # 'alt.pdist': Name of the cumulative distribution function under the alternative distribution. Some choices are
  #    \code{pnorm, pexp, pcauchy, plaplace, pt, pchisq, pf, punif, pbinom, pgeo, ppois.}
  # ...: Optional arguments to \code{alt.pdist}, EXCLUDING the first argument.
  # Examples:
  #      # Alternative distribution is Normal.
  #      power.binom.test(30,0.05,"greater",55,pnorm,55.7,2.5)
  #
  #      # Alternative distribution is Laplace (double exponential).
  #      power.binom.test(30,0.05,"greater",55,plaplace,55.7,2.5)
  if (!is.numeric(n))  stop("'n' must be numeric.") 
  if (n < 1 )  stop("'n' must be a positive integer.") 
  if (!is.numeric(alpha))  stop("'alpha' must be numeric.") 
  if (alpha<=0 | alpha>=1)  stop("'alpha' must be between 0 and 1.") 
  alternative <- abbreviation(alternative, c("two.sided", "less", "greater"))
  if ( !(alternative %in% c("two.sided", "less", "greater")) ) 
       stop("'alternative' must be 'two.sided', 'less', or 'greater'.")
  if (!is.numeric(null.median))  stop("'null.median' must be numeric.") 
  if (length(null.median)!=1)  stop("'null.median' must be scalar.")
  if (!is.function(alt.pdist))  stop("'alt.pdist' must be the alternative distribution.")
  if (alternative=="two.sided") {
     q <- qbinom(alpha/2, n, 0.5) ;    q <- q - (pbinom(q,n,0.5)>alpha/2)
     prob <- pbinom( q, n, 1-match.fun(alt.pdist)(null.median, ...) )  +
             pbinom( q, n, match.fun(alt.pdist)(null.median, ...) ) }
  if (alternative=="less") {
     q <- qbinom(alpha, n, 0.5) ;    q <- q - (pbinom(q,n,0.5)>alpha)
     prob <- pbinom( q, n, 1-match.fun(alt.pdist)(null.median, ...) ) }
  if (alternative=="greater") {
     q <- qbinom(alpha, n, 0.5) ;    q <- q - (pbinom(q,n,0.5)>alpha)
     prob <- pbinom( q, n, match.fun(alt.pdist)(null.median, ...) ) }
  prob
}
