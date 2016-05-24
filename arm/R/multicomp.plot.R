#==============================================================================
# Multiple Comparison Plot
#==============================================================================
multicomp.plot <- function(object, alpha=0.05, main = "Multiple Comparison Plot", 
      label = NULL, shortlabel = NULL, show.pvalue = FALSE,
      label.as.shortlabel = FALSE, label.on.which.axis = 3,
      col.low =  "lightsteelblue", col.same =  "white", col.high = "lightslateblue",
      vertical.line = TRUE, horizontal.line = FALSE,
      vertical.line.lty = 1, horizontal.line.lty = 1, mar=c(3.5,3.5,3.5,3.5)) 
{

  # object check: S4  methods instead?!
  if (!is.data.frame(object)){
    if(is.matrix(object)){
        object <- as.data.frame(object)
    }
    else stop ( message = "object must be a matrix or a data.frame" )
  }
  ind <- dim( object ) [2]
  name <- dimnames( object ) [[2]]
  # label
  if( is.null( label ) ) {
    label <- name
  } else if( length( label ) != ind ) {
    stop( message = "you must specify all the label" )
  }
  # short label
  if( !is.null( shortlabel ) && length( shortlabel ) != ind ){
    stop( message = "you must specify all the short label" )
  } 
  else if( is.null( shortlabel ) && label.as.shortlabel  ){
    shortlabel <- abbreviate( label, minlength = 2)
  } 
  ################################
  # Calculate bayesian p-value
  ################################
  bayes.pvalue <- matrix( 0, ind, ind )
  bayes.signif <- matrix( 0, ind, ind )
  for( i in 1:ind ) {
    for( j in 1:ind ) {
      bayes.pvalue[i, j] <- .pvalue( object[ , j], object[ , i] )
    }
  }
  for( i in 1:ind ) {
    for( j in 1:ind ) {
      bayes.signif[i, j] <- .is.significant( bayes.pvalue[i, j], alpha = alpha ) 
    }
  }
  dimnames( bayes.pvalue ) <- list( label, label )
  diag( bayes.signif ) <- 0
  dimnames( bayes.signif ) <- list( label, label )
  bayes.signif <- bayes.signif [  , ind:1]
  bayes.pvalue <- bayes.pvalue [  , ind:1]
  ################################
  # Plot
  ################################          
  maxchar <- max(sapply(label, nchar))
  mar.idx <- label.on.which.axis
  
  par(mar=mar)
  min.mar <- par('mar')
  if(mar.idx==3){
    mar[mar.idx] <- min(min.mar[mar.idx], trunc(mar[mar.idx] + maxchar/3)) + mar[mar.idx] + 0.1
  }
  else {
    mar[mar.idx] <- min(min.mar[mar.idx], trunc(mar[mar.idx] + maxchar/2)) + 0.1
  }
  par(mar=mar)
  image( 1:nrow( bayes.signif ), 1:ncol( bayes.signif ), 
          bayes.signif, ylab = "", xlab = "", yaxt = "n", xaxt = "n",
          col = c( col.low, col.same, col.high ) )
  box( "plot" )
  axis(2, at = 0, labels = "", las = 1, line = 0, tick = FALSE, 
          xaxs = "i", yaxs = "i" )
  axis(mar.idx, at = 1:nrow( bayes.signif ),line = -0.8, las = 2 , cex = 0.3,
          labels = label, tick = FALSE, xaxs = "i") 
  title( main = main, line = mar[3] - 3 )
  
  for( a in 1:ind ) {
    if( vertical.line ) {
      lines( c( a + 0.5, a + 0.5 ), c( 0, ind + 1 ), lty =  vertical.line.lty )
    }
    if( horizontal.line ) {
      lines(  c( 0, ind + 1 ), c( a + 0.5, a + 0.5 ), lty = horizontal.line.lty )
    }
    if( !is.null( shortlabel ) ) {
      for( b in 1:ind ) {
        if( show.pvalue ){
          text( a, b, ( round( bayes.pvalue, 2 ) )[a,b], cex = 0.5 )
        } else {
          text(  a, b, shortlabel[ind+1-b], cex = 0.7 )
        }
      }
    }
  }
  invisible( list( pvalue = bayes.pvalue, significant = bayes.signif ) )
}

mcplot <- multicomp.plot
