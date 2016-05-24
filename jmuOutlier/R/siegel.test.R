siegel.test <-
function(x, y, alternative=c("two.sided", "less", "greater"), 
                        reverse=FALSE, all.perms=TRUE, num.sim=20000)   {
  # Performs the Siegel-Tukey test, where ties are handled by averaging ranks,
  #   not by asymptotic approximations.
  # 'x': Numeric vector of data values.
  # 'y': Numeric vector of data values.
  # 'alternative': A character string specifying the alternative hypothesis, and
  #   must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}.  
  #   Only the initial letter needs to be specified.
  # 'reverse': Logical; If \code{FALSE} (default), then assign rank 1 to the smallest observation.
  #                     If \code{TRUE}, then assign rank 1 to the largest observation.
  # 'all.perms': Logical.  The exact p-value is attempted when \code{all.perms} (i.e., all permutations) 
  #   is \code{TRUE} (default), and is simulated when \code{all.perms} is \code{FALSE} or when
  #   computing an exact p-value requires more than \code{num.sim} calculations.
  # 'num.sim': The upper limit on the number of permutations generated.
  # Example:  siegel.test( c(13, 34, 2, 19, 49, 63), c(17, 29, 22) )
  # Example:  siegel.test( c(13, 34, 2, 19, 49, 63), c(17, 29, 22), reverse=TRUE )
  if ( !is.numeric(x) )   stop("'x' must be numeric.")
  if ( !is.numeric(y) )   stop("'y' must be numeric.")
  alternative <- abbreviation(alternative, c("two.sided", "less", "greater"))
  if ( !(alternative %in% c("two.sided", "less", "greater")) ) 
       stop("'alternative' must be 'two.sided', 'less', or 'greater'.")
  if ( !is.logical(reverse) )   stop("'reverse' must be logical.")
  if ( !is.logical(all.perms) )   stop("'all.perms' must be logical.")
  if ( !is.numeric(num.sim) )  stop("'num.sim' must be numeric.")
  num.sim = floor(num.sim) ;    if ( num.sim < 1 )  stop("'num.sim' must be at least 1.")
  z <- sort(c(x,y)) ;     lz <- length(z)
  lower <- 1;  upper <- lz;  j <- rep(NA, lz)
  for (i in 1:lz)   {
    if (i-4*as.integer(i/4)<=1)  { j[lower] <- i ;   lower <- lower+1 }
    else {j[upper] <- i;   upper <- upper-1}      }
  if (reverse) { j <- rev(j) }
  z0 <- z[order(j)]
  rank.x.mat <- matrix(NA, length(x), lz)
  rank.y.mat <- matrix(NA, length(y), lz)
  for (k in 1:lz)  {
    for (i in 1:length(x))  { if (x[i]==z0[k]) {rank.x.mat[i,k] <- k} }
    for (i in 1:length(y))  { if (y[i]==z0[k]) {rank.y.mat[i,k] <- k} }    }
  rank.x <- apply(rank.x.mat,1,mean,na.rm=TRUE)
  rank.y <- apply(rank.y.mat,1,mean,na.rm=TRUE)
  if (alternative=="two.sided")     
    {p.value <- perm.test(rank.x,rank.y,alternative="two.sided",all.perms=all.perms,
                          num.sim=num.sim,stat=mean)$p.value}
  if (alternative=="less")
    {p.value <- perm.test(rank.x,rank.y,alternative="greater",all.perms=all.perms,
                          num.sim=num.sim,stat=mean)$p.value}
  if (alternative=="greater")
    {p.value <- perm.test(rank.x,rank.y,alternative="less",all.perms=all.perms,
                          num.sim=num.sim,stat=mean)$p.value}
  structure(list("   Siegel-Tukey test", alternative=alternative,
                  rank.x=rank.x, rank.y=rank.y, p.value=p.value))
}
