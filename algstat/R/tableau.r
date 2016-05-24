#' Tableau Notation for Markov
#'
#' Print the tableau notation for a Markov move.  See the reference provided, p. 13.
#' 
#' @param move a markov move matrix, where the columns are moves in vector form (e.g. the output of markov)
#' @param dim the dimensions of the table form of the move, oftentimes a vector of the number of levels of each variable in order
#' @return an object of class tableau
#' @export tableau
#' @references Drton, M., B. Sturmfels, and S. Sullivant (2009). \emph{Lectures on Algebraic Statistics}, Basel: Birkhauser Verlag AG.
#' @examples
#' \dontrun{
#' 
#' # 2x2 independence example
#' # following convention, the first index indicates rows
#' varlvls <- c(2,2)
#' facets <- list(1,2)
#' ( A <- hmat(varlvls, facets) )
#' markov(A)
#' markov(A, "vec")
#' markov(A, "tab", varlvls)
#' markov(A, "tab", varlvls, TRUE)
#' tableau(markov(A), varlvls)
#' 
#' 
#' 
#' 
#' 
#' 
#' # LAS example 1.2.12, p.17  (no 3-way interaction)
#' varlvls <- c(2,2,2)
#' facets <- list(c(1,2), c(1,3), c(2,3))
#' ( A <- hmat(varlvls, facets) )
#' markov(A)
#  tableau(markov(A), varlvls)
#' 
#'
#' 
#' 
#'
#' }
#' 
tableau <- function(move, dim){

  if(is.vector(move) && (is.integer(move) || is.numeric(move))) 
	  move <- t(t(move))

  if(ncol(move) == 1){  
    p <- length(dim)
    
    move <- t(t(tab2vec(vec2tab(move, dim))))
    row.names(move) <- str_replace_all(row.names(move), ",", "")
	
    t1 <- row.names(move)[move > 0]
    t1 <- matrix(
      unlist(lapply(strsplit(t1, ''), as.numeric)),
      ncol = p, byrow = TRUE
    )
  
    t2 <- row.names(move)[move < 0]
    t2 <- matrix(
      unlist(lapply(strsplit(t2, ''), as.numeric)),
      ncol = p, byrow = TRUE
    )
  
    out <- list(t1, t2) 
    
    
  } else {
  	
  	
    out <- list()
    for(k in 1:ncol(move)){
      out[[k]] <- tableau(move[,k,drop = FALSE], dim = dim)  
    }
    
  }
  
  
  
  class(out) <- 'tableau'
  out
}














#' Pretty printing of tableau output.
#'
#' Pretty printing of tableau output.
#' 
#' @param x an object of class tableau
#' @param ... ...
#' @usage \method{print}{tableau}(x, ...)
#' @return Invisible string of the printed object.
#' @export
#' @examples
#' 
#' # see ?tableau
#' 
#' 
print.tableau <- function(x, ...){
  if((length(x) == 2) && all(sapply(x, is.matrix))){
    plus <- x[[1]]
    minus <- x[[2]]
    p <- ncol(x[[1]])
    changes <- nrow(x[[1]])  
    for(k in 1:changes){
      if(k == floor(changes/2)){
        sep <- ' - '
      } else {
        sep <- '   '
      }
      line2print <- paste(
        paste(plus[k,], collapse = ' '),
        sep,
        paste(minus[k,], collapse = ' '),'\n'
      )
      cat(line2print)
    }  	
  } else {
    for(k in 1:length(x)){
      print.tableau(x[[k]])
      if(k < length(x)) cat('\n')
    }
  }
}
