#' Check if R is in the startup sequence.
#' 
#' @export
is_r_startup<- function(){
    root <- sys.call(1)
    !is.null(root) && (deparse(root) == ".First.sys()")
}
