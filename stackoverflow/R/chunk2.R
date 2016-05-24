#' Split a vector into n chunks
#' 
#' @param x a vector
#' @param n number of chunks
#' 
#' @author 
#'   \href{http://stackoverflow.com/users/1563634/mathheadinclouds}{mathheadinclouds}, 
#'   \href{http://stackoverflow.com/users/1737569/dis-shishkov}{Dis Shishkov}
#' @references \url{http://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r}
#' @export
#' @examples
#'  chunk2(1:30, 6)

chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

