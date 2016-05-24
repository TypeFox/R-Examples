perm.cor.test <-
function(x, y=NULL, alternative=c("two.sided", "less", "greater"), 
                          method=c("pearson", "spearman"), num.sim=20000 ) {
  # A permutation correlation test is performed.
  # P-values are approximated based on randomly sampling the permutations,
  #   where \code{x} and \code{y} are the paired vectors of observations.
  # 'x': Numeric vector of design variable.
  # 'y': Numeric vector of response variable.
  # 'alternative': A character string specifying the alternative hypothesis, and
  #   must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}.  
  #   Only the initial letter needs to be specified.
  # 'num.sim': The upper limit on the number of permutations generated.
  # `method' is a character string specifying the type of correlation, and
  #   must be one of \code{"pearson"} (default) or \code{"spearman"}.  
  #   Only the initial letter needs to be specified.
  # `num.sim' is the number of permutations generated.
  # Examples:   x = c( 4, 6, 8, 11 ) ;   y = c( 19, 44, 15, 13 )
  #             perm.cor.test( x, y, "less", "pearson" ) 
  #             perm.cor.test( x, y, "less", "spearman" ) 
  if ( !is.numeric(x) )   stop("'x' must be numeric.")
  if ( !is.numeric(y) & !is.null(y) )   stop("'y' must be numeric or NULL.")
  alternative <- abbreviation(alternative, c("two.sided", "less", "greater"))
  if ( !(alternative %in% c("two.sided", "less", "greater")) ) 
       stop("'alternative' must be 'two.sided', 'less', or 'greater'.")
  method <- abbreviation(method, c("pearson", "spearman"))
  if ( !(method %in% c("pearson", "spearman")) ) 
       stop("'method' must be 'pearson' or 'spearman'.")
  if ( !is.numeric(num.sim) )  stop("'num.sim' must be numeric.")
  num.sim = floor(num.sim);   if ( num.sim < 1 )  stop("'num.sim' must be at least 1.")
  if (is.null(y))   {y <- x[,2] ;  x <- as.numeric(x[,1])}
  if (length(x)!=length(y)) stop("'x' and 'y' must have the same length.")
  if (method=="spearman")  {x <- rank(x);  y <- rank(y)}
  if (method %in% c("pearson", "spearman")) {
    test.stat0 <- cor(x,y);  test.stat <- NULL
    for (i in 1:num.sim)  {test.stat <- c(test.stat, cor(x, sample(y)))}
    if (alternative=="two.sided") p.value <- mean(abs(test.stat)>=abs(test.stat0))
    if (alternative=="less") p.value <- mean(test.stat<=test.stat0)
    if (alternative=="greater") p.value <- mean(test.stat>=test.stat0)   }
  output1 <- paste("Permutation correlation test.  Method is", method)
  output2 <- paste("p-value was estimated based on", num.sim, "simulations.")
  structure(list(output1, output2, alternative=alternative, p.value=p.value))
}
