#####
#
#  Check what is ready

whatIsSpecified <- function(data){

  N <- length(data)  
  res <- list()
  for(i in 1:N){
    res[[i]] <- list()
    res[[i]]$x <- TRUE
    res[[i]]$y <- TRUE
    res[[i]]$lambda <- TRUE
    res[[i]]$sigma <- TRUE
    res[[i]]$SB <- TRUE
    res[[i]]$smoothed <- TRUE
    
    if( is.null(data[[i]]$x) || any(is.na(data[[i]]$x)) ||  !any(data[[i]]$x!=0) )
      res[[i]]$x <- FALSE  
    if( is.null(data[[i]]$y) || any(is.na(data[[i]]$y)) ||  !any(data[[i]]$y!=0) )
      res[[i]]$y <- FALSE  
    if( is.null(data[[i]]$lambda) || any(is.na(data[[i]]$lambda)) ||  !any(data[[i]]$lambda!=0) )
      res[[i]]$lambda <- FALSE  
    if( is.null(data[[i]]$sigma) || any(is.na(data[[i]]$sigma)) ||  !any(data[[i]]$sigma!=0) )
      res[[i]]$sigma <- FALSE  
    if( is.null(data[[i]]$SB) || any(is.na(data[[i]]$SB)) ||  !any(data[[i]]$SB!=0) )
      res[[i]]$SB <- FALSE  
    if( is.null(data[[i]]$smoothed) || any(is.na(data[[i]]$smoothed)) ||  !any(data[[i]]$smoothed!=0) )
      res[[i]]$smoothed <- FALSE        
  }  
   
  return(res) 
}