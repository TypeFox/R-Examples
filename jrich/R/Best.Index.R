#'
#' @title Jack-knife indices in a single topology m times and evaluates a 
#' success rule.
#'
#' @description The function jack-knifes the terminals and
#' calculates the indices value m (=replicates) times.  
#' 
#' @return The function returns the success that correspond to  
#' obtain the same ranking for X,Y positions, established as the vector
#' success (by default success)).
#'
#' @param tree is a single tree with n terminals, an ape phylo object.
#' 
#' @param distribution species distributions in n areas, a data.frame
#' 
#' @param jtip is the proportion of terminals to delete, real (range 0-1).
#' 
#' @param replicates is the number of replicates, an integer.
#' 
#' @param success the measure of the success, a vector.
#' 
#' @return The function returns the success that corresponds to  
#' obtain the same ranking for X,Y positions, established as the vector
#' success (by default success))
#' 
#' @examples
#' library(jrich)
#' data(tree)
#' data(distribution)
#'
#' Best.Index(tree = tree, distribution = distribution, jtip =0.5, replicates =10, success=1)
#'
#'
#'@author Miranda-Esquivel Daniel R.
#'
#'




Best.Index <-
function (tree = tree, distribution = distribution, jtip = jtip,
          replicates=replicates, success=c(success) ) {
          
  rank <- Rank.Indices(Calculate.Index(tree = tree,distribution = distribution))
  
  aciertos <- NULL
  
  aciertos$I <- aciertos$Ie <- aciertos$Is <- aciertos$Ise <- aciertos$W <- aciertos$We <- aciertos$Ws <- aciertos$Wse <-0
    
  for (i in 1:replicates){
    
    jack <- Rank.Indices(Calculate.Index(tree = tree, distribution = distribution, jtip))
    
    if(all(rank$I[success] == jack$I[success])){
		ok = 1}else{
		ok=0
    }
    
    aciertos$I <- aciertos$I+ok
    
    if(all(rank$Ie[success] == jack$Ie[success])){
		ok = 1}else{
		ok=0
	}
    aciertos$Ie <- aciertos$Ie+ok
    
    if(all(rank$Is[success] == jack$Is[success])){
		ok = 1}else{
		ok=0
	}
    aciertos$Is <- aciertos$Is+ok
    
    if(all(rank$Ise[success] == jack$Ise[success])){
		ok = 1}else{
		ok=0
	}
    aciertos$Ise <- aciertos$Ise+ok
      
    if(all(rank$W[success] == jack$W[success])){
		ok = 1}else{
		ok=0
	}
    aciertos$W <- aciertos$I+ok
    
    if(all(rank$We[success] == jack$We[success])){
		ok = 1}else{
		ok=0
	}
    aciertos$We <- aciertos$We+ok
    
    if(all(rank$Ws[success] == jack$Ws[success])){
		ok = 1}else{
		ok=0
	}
    aciertos$Ws <- aciertos$Ws+ok
    
    if(all(rank$Wse[success] == jack$Wse[success])){
		ok = 1}else{
		ok=0
	}
    aciertos$Wse <- aciertos$Wse+ok      
  }
  
  aciertos <- as.data.frame(aciertos)
  
  aciertos <- aciertos/replicates*100

  return(aciertos)
  
}
