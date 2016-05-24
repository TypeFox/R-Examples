#' Generate common order constraints
#'
#' @description Automatically generates the constraints in the format used by \code{\link{clme}}. Allowed orders are simple, simple tree, and umbrella orders.
#'
#'
#' @param P1 the length of \eqn{\theta_1}{theta_1}, the vector constrained coefficients.
#' @param constraints List with the elements \code{order}, \code{node}, and \code{decreasing}. See Details for further information.
#'
#' @details 
#' The elements of \code{constraints} are:
#' \itemize{
#' \item \code{order}: string. Currently \dQuote{simple}, \dQuote{simple.tree} and \dQuote{umbrella} are supported.
#' \item \code{node}: numeric, the node of the coefficients (unnecessary for simple orders).
#' \item \code{decreasing}: logical. For simple orders, is the trend decreasing? For umbrella and simple tree, does the nodal parameter have the greatest value (e.g., the peak, instead of the valley)?
#' }
#' 
#' See \code{\link{clme}} for more information and a depiction of these three elements.
#' 
#' @return
#' The function returns a list containing the elements of input argument \code{constraints} as well as
#' \itemize{
#' \item{ \code{A} }{matrix of dimension \eqn{r \times 2}{r x 2} containing the order constraints, where r is the number of linear constraints.}
#' \item{ \code{B} }{matrix containing the contrasts necessary for computation of the Williams' type test statistic (may be identical to \code{A}).}
#' \item{ \code{Anull} }{matrix similar to \code{A} which defines all possible constraints. Used to obtain parameter estimates under the null hypothesis.}
#' \item{ \code{order} }{the input argument for \code{constraints\$order}.}
#' \item{ \code{node} }{the input argument for \code{constraints\$node}.}
#' \item{ \code{decreasing} }{ the input argument for \code{constraints\$decreasing}}
#' }
#' See \code{\link{w.stat}} for more information on \code{B}
#' 
#' 
#' @note
#' The function \code{\link{clme}} also utilizes the argument \code{constraints}. For \code{clme}, this argument may either be identical to the argument of this function, or may be the output of \code{create.constraints} (that is, a list containing appropriate matrices \code{A}, \code{Anull}, and if necessary, \code{B}).
#' 
#' An example the the \code{A} matrix might be:
#'   \tabular{ccc}{
#'     [1,] \tab [,1] \tab [,2] \cr
#'     [2,] \tab 1    \tab 2    \cr
#'     [3,] \tab 2    \tab 3    \cr
#'     [4,] \tab 4    \tab 3    \cr
#'     [5,] \tab 5    \tab 4    \cr
#'     [6,] \tab 6    \tab 5    \cr
#'   }
#' This matrix defines what \pkg{CLME} describes as a decreasing umbrella order. The first row defines the constraint that \eqn{\theta_1 \leq \theta_2}{theta_1 <= theta_2}, the second row defined the constraint \eqn{\theta_2 \leq \theta_3}{theta_2 <= theta_3}, the third row defines \eqn{\theta_4 \leq \theta_3}{theta_4 <= theta_3}, and so on. The values are indexes, and the left column is the index of the parameter constrained to be smaller.
#' 
#' 
#' @seealso
#' \code{\link{clme}},
#' \code{\link{w.stat}}
#' 
#' @examples
#' \dontrun{
#'   # For simple order, the node does not matter
#'   create.constraints( P1 = 5, constraints = list( order='simple' , 
#'                                                   decreasing=FALSE ))
#'   
#'   # Compare constraints against decreasing=TRUE
#'   create.constraints( P1 = 5, constraints=list( order='simple' , 
#'                                                 decreasing=TRUE ))
#'   
#'   # Umbrella order
#'   create.constraints( P1 = 5, constraints=list( order='umbrella' , node=3
#'                                                 , decreasing=FALSE ))
#' }
#' 
#' 
#' @importFrom utils combn
#' @export
#' 
create.constraints <- function( P1, constraints ){
  
  Q1 <- P1-1

  
  order      <- tolower(constraints$order)
  node       <- constraints$node
  decreasing <- constraints$decreasing
  
  if( is.null(node) ){
    node <- 1
  }
  
  if( (order %in% c("simple", "simple.tree", "umbrella"))==FALSE ){
    stop("'order' must be one or more of: simple, simple.tree, umbrella")
  }
  
  ## Revert to simple order if umbrella has node at extreme
  if( order=="umbrella" & node %in% c(1,P1) ){
    order <- "simple"
    node  <- 1
    if( node==P1 ){
      if( decreasing==TRUE ){
        decreasing <- FALSE
      } else {
        decreasing <- TRUE
      }
    }
  }
  
  
  A <- matrix( 0, nrow=Q1, ncol=2 )
  
  ## Simple order
  ## e.g. mu_1 <= mu_2 <= ... <= mu_K
  if( order=="simple" ){
    if( decreasing==TRUE ){
      A <- as.matrix(cbind( 1:Q1+1 , 1:Q1 ))
      B <- matrix( c(P1,1) , nrow=1 )
    } else{
      A <- as.matrix(cbind( 1:Q1 , 1:Q1 + 1 ))
      B <- matrix( c(1,P1) , nrow=1 )
    }
    node    <-  NULL
  }
  
  ## Simple tree order
  ## e.g. mu_1 <= mu_i ; i=2,...,P1
  if( order=="simple.tree" ){
    if( decreasing==TRUE ){
      A <- as.matrix( cbind( (1:P1)[-node] , rep(node,Q1) ) )
    } else{
      A <- as.matrix(cbind( rep(node,Q1)  , (1:P1)[-node] ))
    }
    B <- A    
  }
  
  ## Umbrella order
  ## e.g. mu_1 <= mu_2 <= ... <= mu_b >= mu_{b+1} >= ... >= mu_K
  if( order=="umbrella" ){    
    if( decreasing==TRUE ){
      for( ii in 1:(node-1) ){
        A[ii,] <- c(ii,ii+1)
      }
      for( ii in node:Q1 ){
        A[ii,] <- c(ii+1,ii)
      }
      B <- as.matrix( rbind( c(1,node) , c(P1,node) ) )
      
    } else{
      
      for( ii in 1:(node-1) ){
        A[ii,] <- c(ii+1,ii)
      }
      for( ii in node:Q1 ){
        A[ii,] <- c(ii,ii+1)
      }
      B <- as.matrix( rbind( c(node,1) , c(node,P1) ) )      
    }
  }
  
  
  ## Make the null A-matrix
  Anull <- t(combn( 1:P1 , m=2 ))
  Anull <- rbind( Anull, Anull[,2:1] )
  
  ## If A has only one row, activeSet() from package "isotone" causes error.
  ## Placing same constraint twice alleviates this problem.
  if( nrow(A)==1 ){
    A     <- rbind( A    , A     )
    Anull <- rbind( Anull, Anull )
  }
  
  # Return the constraints object
  new_constraints <- list( A = A, B = B, Anull = Anull, order=order, node=node, decreasing=decreasing )
  return(new_constraints)
}
