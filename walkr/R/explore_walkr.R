#' Explore Walkr
#' 
#' This function takes in a list of chains and diagnoses them using
#' the shinyStan interface. The app contains the confidence
#' interval of each dimension's coordiates across the whole sets of points,
#' Gelman-Rubin statistics, trace-plots, and other diagnostic tools for
#' examining convergence.
#' 
#' @param x is the set of points sampled, with its columns as the sampled points.
#'        If multiple chains are present, then the columns should be ordered
#'        such that each chain follow each other. 
#' 
#' @return a shiny interface that display the diagnostics of the MCMC random walk
#' @export
#' @importFrom shinystan launch_shinystan as.shinystan
#' @import ggplot2

explore_walkr <- function(x) {
  
  ## checking
  stopifnot(is.list(x))
  
  ## taking the first chain
  this.df <- x[[1]] 
  chains <- length(x)
  
  ## if there are more than 1 chains, we 
  ## make all the chains into an ordered 
  ## matrix to be manipulated below
  
  if(chains > 1){
    
    for (i in 2:chains){
      this.df <- cbind(this.df, x[[i]])
    }
    
  }  
  
  ## iterations is the number of points in hitandrun
  ## parameters is the dimension of the sample space
  
  iterations <- ncol(this.df) / chains
  parameters <- nrow(this.df)
  
  ## we have to transpose here in order to make the dimensions fit 
  
  this.df <- t(this.df)
  
  ## remake the data set into 3 dimensional array
  
  dim(this.df) <- c(iterations, chains, parameters)
  
  ## using the stan package
  
  shinystan::launch_shinystan(shinystan::as.shinystan(this.df))
  
}