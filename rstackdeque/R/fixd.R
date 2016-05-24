#' @export
#' @title Fix an rdeque
#' 
#' @description Maintains the invariant that there is always something in two stacks used by 
#' rdeques under the hood so long as there is 2 more elements in the rdeque.
#' 
#' @details In fact, fix will be called whenever there are fewer than 6 elements in both
#' the front and end of the deque. Generally this method is \eqn{O(N)}, and so a full copy is returned.
#' @param d rdeque to fix.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @references Okasaki, Chris. Purely Functional Data Structures. Cambridge University Press, 1999.
#' @return fixed, "balanced" deque.
fixd <- function(d, ...) {UseMethod("fixd", d)}