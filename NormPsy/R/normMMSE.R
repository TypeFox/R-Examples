normMMSE <- function(x)
{
 if(missing(x)) stop("Please specify the MMSE score to transforme in argument 'x'")  
 if(!all(na.omit(x) %in% 0:30)) stop("The score to transforme must be an integer between 0 and 30")  
  
 transfY <- c(0, 2.91, 5.48, 7.76, 9.77, 11.57, 13.19, 14.67, 16.05, 17.37, 
  18.68, 20.01, 21.38, 22.83, 24.39, 26.07, 27.91, 29.93, 32.17, 
  34.64, 37.37, 40.4, 43.74, 47.4, 51.44, 55.98, 61.18, 67.25, 
  74.61, 84.32, 100) 
  
  return(transfY[x+1])
}

