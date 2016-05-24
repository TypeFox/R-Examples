rmd.test <-
function(x, y, alternative=c("two.sided", "less", "greater"), all.perms=TRUE,
                     num.sim=20000) {
  # A permutation test is performed based on the estimated rmd,
  #   the ratio of the absolute mean deviances.
  # 'x': Numeric vector of data values.
  # 'y': Numeric vector of data values.
  # 'alternative': A character string specifying the alternative hypothesis, and
  #   must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}.  
  #   Only the initial letter needs to be specified.
  # 'all.perms': Logical.  The exact p-value is attempted when \code{all.perms} (i.e., all permutations) 
  #   is \code{TRUE} (default), and is simulated when \code{all.perms} is \code{FALSE} or when
  #   computing an exact p-value requires more than \code{num.sim} calculations.
  # 'num.sim': The upper limit on the number of permutations generated.
  # Example:  rmd.test( c(13, 34, 2, 19, 49, 63), c(17, 29, 22) )
  # Example:  rmd.test( c(13, 34, 2, 19, 49, 63), c(17, 29, 22), "greater" )
  if ( !is.numeric(x) )   stop("'x' must be numeric.")
  if ( !is.numeric(y) )   stop("'y' must be numeric.")
  alternative <- abbreviation(alternative, c("two.sided", "less", "greater"))
  if ( !(alternative %in% c("two.sided", "less", "greater")) ) 
       stop("'alternative' must be 'two.sided', 'less', or 'greater'.")
  if ( !is.logical(all.perms) )   stop("'all.perms' must be logical.")
  if ( !is.numeric(num.sim) )  stop("'num.sim' must be numeric.")
  num.sim = floor(num.sim);   if ( num.sim < 1 )  stop("'num.sim' must be at least 1.")
  max.loop <- 10;        l1 <- min(length(x), length(y)) ;    tol <- 1e-10
  alternative <- abbreviation(alternative, c("two.sided", "less", "greater"))
  l2 <- max(length(x), length(y));    l3 <- l1+l2
  x <- x-median(x) ;   y <- y-median(y)
  rmd <- function(x, y, alternative) {
    test <- mean(abs(x)) / mean(abs(y))
    if (alternative=="two.sided") {test <- max(test, 1/test)}
    return( test )      }
  test.stat0 <- rmd(x, y, alternative) ;    test.stat <- NULL
  if (choose(l3,l1)>num.sim) all.perms <- FALSE
  if (all.perms & l1<=max.loop) {
    output <- "Exact p-value was calculated."
    index <- rep(NA, max(l3, max.loop))
    for (i1 in 1:(l2+1))  {
    for (i2 in (i1+1):max( i1+1, (l2+2)*(l1>=2) ))  {
    for (i3 in (i2+1):max( i2+1, (l2+3)*(l1>=3) ))  {
    for (i4 in (i3+1):max( i3+1, (l2+4)*(l1>=4) ))  {
    for (i5 in (i4+1):max( i4+1, (l2+5)*(l1>=5) ))  {
    for (i6 in (i5+1):max( i5+1, (l2+6)*(l1>=6) ))  {
    for (i7 in (i6+1):max( i6+1, (l2+7)*(l1>=7) ))  {
    for (i8 in (i7+1):max( i7+1, (l2+8)*(l1>=8) ))  {
    for (i9 in (i8+1):max( i8+1, (l2+9)*(l1>=9) ))  {
    for (i10 in (i9+1):max( i9+1, (l2+10)*(l1>=10) ))  {
      index <- c(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10)
      index <- c( index[1:l1], setdiff(c(1:l3),  index[1:l1]) )
      test.stat <- c(test.stat, rmd( c(x,y)[index[1:l1]], c(x,y)[index[(l1+1):l3]], alternative ) )
    }}}}}}}}}} # Number of parentheses is max.loop
  } # End of if
  else { # Use simulations.
    output <- paste("p-value was estimated based on", num.sim, "simulations.")
    for (i in 1:num.sim) {
      index <- sample(l3)
      test.stat <- c(test.stat, rmd( c(x,y)[index[1:l1]], c(x,y)[index[(l1+1):l3]], alternative ) )
    }
  }
  if (length(x)>length(y) & alternative!="two.sided") { test.stat <- 1 / test.stat }
  if (alternative=="less") {
    p.value <- mean( test.stat - test.stat0 <= tol ) }
  else {
    p.value <- mean( test.stat - test.stat0 >= -tol) }
  structure(list(output, alternative=alternative, rmd.hat=test.stat0, p.value=p.value))
}
