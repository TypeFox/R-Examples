# ==============================================================================
# residual plot for the observed values
# ==============================================================================
residual.plot <- function ( Expected, Residuals, sigma, 
                            main = deparse(substitute( Expected )), 
                            col.pts = "blue", col.ctr = "red", 
                            col.sgm = "black", cex = 0.5, gray.scale = FALSE, 
                            xlab="Predicted", ylab="Residuals", ... ) {
  if( gray.scale == TRUE ) { 
    col.pts <- "black";
    col.ctr <- "black";
    col.sgm <- "gray60";
  }
  plot( Expected[!is.na( Residuals )], Residuals[ !is.na( Residuals ) ],
         xlab = xlab, ylab = ylab, main = main, col = col.pts,
          pch = 19, cex = cex, ... );
  #mtext( "Residuals vs Predicted", 3, cex= 0.6 )  #, adj=1 );
  # add the zero line for clarity
  abline ( h = 0, lty = "dashed", col = col.ctr );
  # residual s.e.
  resid.se <- sigma;
  # Add two-standard-error lines
  abline ( h =  2*resid.se, lty = "dashed", col = col.sgm );
  abline ( h = -2*resid.se, lty = "dashed", col = col.sgm );
} 
