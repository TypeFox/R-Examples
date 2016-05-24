lineGraph <-
function(x, freq=TRUE, prob=NULL, col="red", ...) {
  # Constructs a line graph.
  # 'x': Vector of numerical observations to be graphed.
  # 'freq': Logical.  If \code{freq} is \code{TRUE} (default) and \code{prob} does not sum to 1, then frequencies (counts) are plotted.
  # If \code{freq} is \code{FALSE} or \code{prob} sums to 1, then relative frequencies are plotted.
  # 'prob': Vector of the probabilities or weights on \code{x},
  #   and does not need to sum to one.  If \code{prob} is \code{NULL}, 
  #   then all \code{x} values are equally weighted.
  # 'col': The color of the plotted lines.  Type \code{colors()} for selections.
  # ...: Optional arguments to \code{\link[graphics]{plot}}.
  # example:  lineGraph( c( rep(6,4), rep(9,7), rep(3,5), 5, 8, 8 ) )
  # example:  lineGraph( c( rep(6,4), rep(9,7), rep(3,5), 5, 8, 8 ), freq=FALSE )
  # example:  lineGraph( 11:14, , c( 12, 9, 17, 5 ) )
  if (!is.numeric(x))  stop("'x' must be numeric.")
  x.old <- x;  x <- union(x.old,x.old);  
  if (is.null(prob)) { for (i in 1:length(x)) prob <- c( prob, sum(x.old==x[i]) ) }
  if (!freq) prob <- prob/sum(prob)
  if (sum(prob)==1 || !freq) {plot( x, prob/sum(prob), ylim=c(0,max(prob)/sum(prob)), 
     type="h", ylab="RELATIVE  FREQUENCY", col=col, ... )}
  else {plot( x, prob, ylim=c(0,max(prob)), type="h", ylab="FREQUENCY", col=col, ... )}
}
