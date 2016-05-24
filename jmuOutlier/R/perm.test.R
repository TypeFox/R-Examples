perm.test <-
function(x, y=NULL, alternative=c("two.sided", "less", "greater"), mu=0,
                      paired=FALSE, all.perms=TRUE, num.sim=20000, plot=FALSE, stat=mean, ...) {
  # Performs one-sample and two-sample permutation tests on vectors of data.
  # The one-sample permutation test is based on stat(x-mu);
  #   the two-sample paired permutation test is based on stat(x-y-mu);
  #   and the two-sample unpaired permutation test is based on stat(x-mu)-stat(y).
  # 'x': A (non-empty) numeric vector of data values.
  # 'y': An optional numeric vector data values.
  # 'alternative': A character string specifying the alternative hypothesis, and
  #   must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}.  
  #   Only the initial letter needs to be specified.
  # 'mu': A number indicating the null value of the location parameter 
  #   (or the difference in location parameters if performing a two-sample test).
  # 'paired': Logical, indicating whether or not a two-sample test should be paired,
  #   and is ignored for a one-sample test.
  # 'all.perms': Logical.  The exact p-value is attempted when \code{all.perms} (i.e., all permutations) 
  #   is \code{TRUE} (default), and is simulated when \code{all.perms} is \code{FALSE} or when
  #   computing an exact p-value requires more than \code{num.sim} calculations.
  # 'num.sim': The upper limit on the number of permutations generated.
  # 'plot': Logical. If \code{TRUE}, then plot the histogram of the permutation distribution;
  #   otherwise, list the p-value.  
  #   Setting \code{plot} to code{TRUE} works best when the permutation distribution
  #   is close to continuous.
  # 'stat': Function, naming the test statistic, such as \code{mean} and \code{median}.
  # ... Optional arguments to \code{stat}; 
  #     and is the second argument to \code{stat} when unspecified.
  #     For example, if \code{stat} equals \code{mean}, then the second argument
  #       \code{trim} denotes the fraction (0 to 0.5) of observations to be trimmed 
  #       from each end of \code{x} and \code{y} before the mean is computed. 
  # Examples: 
  #           # One-sample test
  #             x = rnorm(10,0.5) ;  list(x) ;  perm.test( x, stat=median )
  #
  #           # Two-sample unpaired test
  #           y = rnorm(13, 1 ) ;  list(y) ;  perm.test( x, y )
  if ( !is.numeric(x) )  stop("'x' must be numeric.")
  if ( !is.numeric(y) & !is.null(y) )  stop("'y' must be numeric or NULL.")
  alternative <- abbreviation(alternative, c("two.sided", "less", "greater"))
  if ( !(alternative %in% c("two.sided", "less", "greater")) ) 
       stop("'alternative' must be 'two.sided', 'less', or 'greater'.")
  if ( !is.numeric(mu) )  stop("'mu' must be numeric.")
  if ( !is.logical(paired) )  stop("'paired' must be logical.")
  if ( !is.logical(all.perms) )  stop("'all.perms' must be logical.")
  if ( !is.numeric(num.sim) )  stop("'num.sim' must be numeric.")
  num.sim = floor(num.sim);  if ( num.sim < 1 )  stop("'num.sim' must be at least 1.")
  if ( !is.logical(plot) )  stop("'plot' must be logical.")
  if ( !is.function(stat) )  stop("'stat' must be a function.")
  if (!is.null(y) & !paired)   {  # Perform unpaired two-sample permutation test.
    l1 <- min(length(x), length(y));   x <- x-mu
    l2 <- max(length(x), length(y));    l3 <- l1+l2
    test.stat0 <- match.fun(stat)(x, ...)-match.fun(stat)(y, ...);  test.stat <- NULL
    if (choose(l3,l1)>num.sim) all.perms <- FALSE
    output1 <- "Unpaired two-sample permutation test was performed."
    if (all.perms) {
      output2 <- "Exact p-value was calculated."
      combos <- combn( 1:l3, l1 )
      for (i in 1:dim(combos)[2])   {
        index <- c( combos[,i], setdiff( c(1:l3), combos[,i] ) )
        test.stat <- c(test.stat,
              match.fun(stat)(c(x,y)[index[1:l1]], ...)
             -match.fun(stat)(c(x,y)[index[(l1+1):l3]], ...) )      }
    } # End of if
    else { # Use simulations.
      output2 <- paste("p-value was estimated based on", num.sim, "simulations.")
      for (i in 1:num.sim) {
        index <- sample(l3)
        test.stat <- c(test.stat,
              match.fun(stat)(c(x,y)[index[1:l1]], ...)
             -match.fun(stat)(c(x,y)[index[(l1+1):l3]], ...) )    }      }
    if (length(x)>length(y)) test.stat <- -test.stat
    if (!plot) {
      if (alternative=="less") p.value <- mean( test.stat <= test.stat0 )
      if (alternative=="greater") p.value <- mean( test.stat >= test.stat0 )
      if (alternative=="two.sided") p.value <- mean( abs(test.stat) >= abs(test.stat0) )
      structure(list(output1, output2, alternative=alternative, mu=mu, p.value=p.value))   }
    else {hist(test.stat, main="Histogram of Permutation Distribution",
                xlab="Value of test statistic", ylab="Frequency")}
    }     # END of `if (!is.null(y) & !paired)'
  else { # Perform a one-sample or a paired two-sample permutation test on mu.
    if (is.null(y)) { output1 <- "One-sample permutation test was performed." }
    else { x <- x-y;  output1 <- "Paired permutation test was performed." }
    x <- x-mu;  lx <- length(x); test.stat0 <- match.fun(stat)(x,...); test.stat <- NULL
    if (2^lx>num.sim) all.perms <- FALSE
    if (all.perms) {
      output2 <- "p-value was calculated based on all permutations."
      plus.minus.matrix <- matrix(1, 1, lx)
         for (i in lx:1)  {
            temp <- plus.minus.matrix ;     temp[,i] <- -1
            plus.minus.matrix <- rbind( plus.minus.matrix, temp )       }
      for (ii in 1:2^lx)   {
        test.stat <- c(test.stat, match.fun(mean)(x*plus.minus.matrix[ii,],...))    }
    } # End of if
    else { # Use simulations.
      output2 <- paste("p-value was estimated based on", num.sim, "simulations.")
      for (i in 1:num.sim) test.stat <- c(test.stat, match.fun(stat)((x*sample(c(-1,1),lx,T)),...)) }
    if (!plot) {
      if (alternative=="less") { p.value <- mean( test.stat <= test.stat0 ) }
      if (alternative=="greater") { p.value <- mean( test.stat >= test.stat0 ) }
      if (alternative=="two.sided") { p.value <- mean( abs(test.stat) >= abs(test.stat0) ) }
      structure(list(output1, output2, alternative=alternative, mu=mu, p.value=p.value))}
    else {hist(test.stat, main="Histogram of Permutation Distribution",
                xlab="Value of test statistic", ylab="Frequency")}
  } # End of `else'
}
