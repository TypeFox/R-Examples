#' Create a model-type vector template
#' 
#' This is a small helper function which creates a vector template  quick and easily for the \code{PPall()} function. Modify this template as you like.
#'
#' This function tries to guess the model which was applied to each item by using the matrix of threshold parameters. It only discriminates between GPCM and 4-PL model, and returns a character vector of length equal to the number of items, that contains \code{"GPCM"} or \code{"4PL"} entries depending on the structure of the thres matrix.
#' 
#'@param thres A numeric matrix which contains the threshold parameter for each item. NA is allowed - in fact expected!
#'
#' @seealso \link{PPall}
#'
#'@export
#'
#'@author Manuel Reif
#'
#'@example ./R/.examples.R
#'@keywords Person Parameters
findmodel <- function(thres)
{

if(any(thres[1,] != 0))
  {
   thres <- rbind(0,thres)
  }  
  
maxsc <- apply(thres,2,function(x)(length(x) - sum(is.na(x)))-1)
model2est <- ifelse(maxsc > 1,"GPCM","4PL")

model2est  
}
