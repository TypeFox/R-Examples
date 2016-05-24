#' @export
#' @rdname timer
toc <- function(){
  if (exists(".time_Jmisc")){
    out = Sys.time() - get(".time_Jmisc")
    print(out)
    return(invisible(out))
  } else {
    warning(".time_Jmisc not found, please run tic() before using toc()")
  }
}
