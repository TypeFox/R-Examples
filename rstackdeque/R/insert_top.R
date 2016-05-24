#' @export
#' @title Insert an element onto the top of an rstack
#' 
#' @description Insert an element onto the top of an rstack.
#' 
#' @details Runs in \eqn{O(1)} time worst-case. Does not semantically modify the original structure (i.e, this
#' function is "pure").
#' @param s rstack to insert onto.
#' @param e element to insert.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return modified version of the stack with new element at top.
#' @seealso \code{\link{rstack}} for information on rstacks, \code{\link{without_top}} for removal of top elements.
#' @references Okasaki, Chris. Purely Functional Data Structures. Cambridge University Press, 1999.
#' @examples
#' s <- rstack()
#' s <- insert_top(s, "a")
#' s <- insert_top(s, "b")
#' print(s)
#' 
#' s2 <- insert_top(s, "c")
#' print(s2)
#' print(s)
insert_top <- function(s, e, ...) { UseMethod("insert_top", s) }

