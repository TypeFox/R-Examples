#' Calculates length of catenary
#'
#' Takes left and right endpoints and catenary parameters and gives length
#' 
#' @param x0 left endpoint
#' @param x1 right endpoint
#' @param c1 catenary parameter
#' @param c2 catenary parameter
#' @return length of catenary
#' @author Jono Tuke <simon.tuke@@adelaide.edu.au>
#' @export
#' @note February 11 2013
#' @keywords internal
#' @examples
#' getCatLength(x0=-1,x1=1,c1=1,c2=2)
getCatLength <- function(x0,x1,c1,c2){
  L <-  c1 * ( sinh( (x1-c2)/c1) - sinh((x0-c2)/c1) )
  return(L)  
}