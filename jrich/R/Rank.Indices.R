#'
#' @title Rank indices. 
#' 
#' @description Renk indices according to the areas' absolute position. If the index value is empty, the function assigns a dummy position "X0X"
#'
#' @param index.Value a table with indices values.
#' 
#' @return a table with the decreasing order of the areas
#' It presents the ties alphabetically 
#' 
#' @examples
#'  ## get the library
#'  library(jrich)
#'  
#'  ## load the data
#'  data(tree) 
#'  data(distribution) 
#'  
#' Rank.Indices(Calculate.Index(tree=tree, distrib = distribution, verbose=FALSE))
#'
#'@author Miranda-Esquivel Daniel R.
#'
#'




Rank.Indices <-
function (index.Value=index.Value) {
 
	ranking <- NULL
	areas <- length(index.Value$area)
 
 
	if (max(index.Value$I) > 0){
		ranking$I <- index.Value$area[order(index.Value$I, decreasing = T)]
		}else{
		ranking$I[1:areas] <- "X0X"  
	}
  
  
	if (max(index.Value$Ie) > 0){
		ranking$Ie <- index.Value$area[order(index.Value$Ie, decreasing = T)]
		}else{
		ranking$Ie[1:areas] <- "X0X"  
	}
 
 
	if (max(index.Value$Is) > 0){
		ranking$Is <- index.Value$area[order(index.Value$Is, decreasing = T)]
		}else{
		ranking$Is[1:areas] <- "X0X"  
	}
 
	if (max(index.Value$Ise) > 0){
		ranking$Ise <- index.Value$area[order(index.Value$Ise, decreasing = T)]
		}else{
		ranking$Ise[1:areas] <- "X0X"  
	}
 
	if (max(index.Value$W) > 0){
		ranking$W <- index.Value$area[order(index.Value$W, decreasing = T)]     
		}else{
		ranking$W[1:areas] <- "X0X"  
	}
 
	if (max(index.Value$We) > 0){
		ranking$We <- index.Value$area[order(index.Value$We, decreasing = T)]  
		}else{
		ranking$We[1:areas] <- "X0X"  
	}
 
	if (max(index.Value$Ws) > 0){
		ranking$Ws <- index.Value$area[order(index.Value$Ws, decreasing = T)]   
		}else{
		ranking$Ws[1:areas] <- "X0X"  
	}
 
	if (max(index.Value$Wse) > 0){
		ranking$Wse <- index.Value$area[order(index.Value$Wse, decreasing = T)]
		}else{
		ranking$Wse[1:areas] <- "X0X"  
	}
 
	ranking$rich <- index.Value$area[order(index.Value$rich, decreasing = T)]
 
	ranking$endem <- index.Value$area[order(index.Value$endem, decreasing = T)]
 
 return(ranking)
 
}
