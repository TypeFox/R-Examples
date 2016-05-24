#' Core IDs
#' 
#' Take the series ID's and convert into a useable format. This produces series ID in a form that can be compared to pith offset file. 
#' @param cores A vector containing core Ids. 
#' @param site.id A vector specifying the position of start and end for the site ID.
#' @param tree.id A vector specifying the position of start and end for the tree ID.
#' @param core.id A vector specifying the position of start and end for the core ID. 
#' @return A vector of core ids, re-formatted. 
#' @export
function.id	<- function( cores, site.id, tree.id, core.id){
  site	<- substr( cores, site.id[ 1 ], site.id[ 2 ])
  tree	<- substr( cores, tree.id[ 1 ], tree.id[ 2 ])
  series	<- substr( cores,core.id[ 1 ], core.id[ 2 ])
  cores.x 	<- toupper(paste("HUP",tree,series,sep="",collapse=NULL))
  return( cores.x )
}

#' Sum non Na's
#' 
#' Count the number of non na values. This is an internal function used in an apply function.  
#' @param x Vector to be counted containing Na's
#' @return The sum of the non na's values. 
#' @export
SumNotNa  <-	function( x ) {
  sum( is.na( x ) == FALSE)
}


#' Exact Confidence Intervals
#' 
#' Calculate the exact confidence intervals for proportions, based in Pearson's approach. 
#' @param p.adjust A vector of concordance indices, these are proportions for which you want to calculate the confidence interval.
#' @param size.x A vector of the number of series used to calculate each index in matrix X. 
#' @param size.y A vector of the number of series used to calculate each index in matrix Y.
#' @param i A recursive variables (i.e. use this function in Apply.)
#' @return  A vector with 2 elements indicating the lower and upper confidence intervals for the concordance indices. 
#' @export
exact.ci<- function(p.adjust, size.x, size.y, i){
  n <- (size.x[i]*size.y[i])
  s <- as.integer(as.vector(p.adjust[i]*n))
  ci.j <- matrix(NA, nrow = length(s), ncol = 2)
  for ( t in 1: length(s)){
    if(is.na(s)){ci.j[,t] <- NA}else{
      if(s[t] > n){s[t] = n}
      ci.j[t,] <- binom.test(s[t],n , p = 1,
                             alternative = c("two.sided"),
                             conf.level = 0.95)$conf.int
    }
  }
  return(ci.j)
}

#' Master ID
#' 
#' Produces series ID in a form that can be compared to pith offset file for the master chronology.
#' @param cores A vector containing core Ids
#' @param site.id A vector specifying the position of start and end for the site ID. 
#' @param tree.id A vector specifying the position of start and end for the tree ID.
#' @param core.id A vector specifying the position of start and end for the core ID.
#' @return A vector of core ids, re-formatted. 
#' @export  
function.id.master	<- function( cores, site.id, tree.id, core.id){
  site	<- substr( cores, site.id[ 1 ], site.id[ 2 ])
  tree	<- substr( cores, tree.id[ 1 ], tree.id[ 2 ])
  series	<- substr( cores,core.id[ 1 ], core.id[ 2 ])
  cores.x 	<- toupper(paste(site,tree,series,sep="",collapse=NULL))
  return( cores.x )
}

