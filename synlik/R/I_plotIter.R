
.plotIter <- function(mat, trans = NULL, type = "line", addPoints = NULL, addPlot = NULL, ...)
{
  if( !(type %in% c("line", "hist")) ) stop("type should be either line or hist")
  
  parNames <- colnames(mat)
  
  if( nrow(mat) > 0 ) 
  {
    nPar <- ncol(mat)
    panelDim <- min( ceiling(sqrt(nPar)), 3 )
    par(mfrow = c(panelDim, panelDim))
    
    # Transform each column if needed
    mat <- .transMatrix(mat, trans)
    
    if( !is.null(addPoints) ){
      xpoints <- addPoints$x
      ypoints <- addPoints$y
      ypoints <- .transMatrix(ypoints, trans)
    }
    
    counter <- 1
    # Plot each column
    for(nam in parNames){
      
      # Either lines
      if(type == "line")
      {
        if( is.null(addPoints) )
        {
        plot(1:nrow(mat), mat[ , nam], type = 'l', main = nam,
             ylab = nam, xlab = "Iteration", ...)
        } else {
          plot(xpoints, ypoints[ , nam], main = nam,
               ylab = nam, xlab = "Iteration",
               ylim =  c(min(c(mat[ , nam], ypoints[ , nam])), max(c(mat[ , nam], ypoints[, nam]))), ...)
          lines(1:nrow(mat), mat[ , nam], col = 2, lwd = 2)
        }
        
      } else {
        # Or histograms
        hist(mat[ , nam],  main = nam, ylab = nam, xlab = nam, ...)
      }
      
      if( !is.null(addPlot) ) get(addPlot)(nam, ...)
      
      if( !(counter %% (panelDim^2)) && (counter != nPar) ) readline(prompt = "Press <Enter> to see the next plot...") 
      counter <- counter + 1
    }
  }
  
  return(NULL)   
}