##' A little stack.
##'
##' TODO: implement as reference class?
##' Note: \code{pop} only removes elements. To retrieve them, use \code{peek}.
##' 
##' @param x the object
##' @param an attribute holding the stack
##' @param n numer of element to peek at and numer of elements to pop (delete), respectively
##' @return \code{push} and \code{pop}: the object with stack in list \code{an} pushed/popped by the \code{n} elements
##'
##' \code{peek}: the \code{n}th stack element (without popping!)
##' @author Claudia Beleites
##' @rdname stack
peek <- function (x, an, n = 1L){
  stack <- attr (x, an)
  if (length (stack) == 0) return (NULL)

  if (length (stack) < n){
    n <- length (stack)
    warning ("n set to ", n)
  }

  stack [[n]]
}

##' @rdname stack
pop <- function (x, an, n = 1L){
  stack <- attr (x, an)
  if (length (stack) == 0) return (x)

  if (length (stack) < n){
    n <- length (stack)
    warning ("n set to ", n)
  }

  if (length (stack) == n)
    attr (x, an) <- NULL
  else
    attr (x, an) <- stack [- seq_len (n)]
      
  x
}

##' @usage push (x, an) <- value
##' @rdname stack
##' @param value list of things to push on the stack.
`push<-` <- function (x, an, value){
  stack <- attr (x, an)
  value <- rev (value)                  # expect: last argument first out
  stack <- c (value, stack)             # prepend the new elements

  attr (x, an) <- stack
      
  x
}

