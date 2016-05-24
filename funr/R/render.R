


## --- showing the output of the returned..

# if error, ultimately, die with a error
# show a clear, clean error message

#' Render output of functions in a print friendly format
#'
#' @description
#' If the function returns with invisible, output is suppressed
#'
#' @param x a output from \code{funr}
#' @param max_rows In case output of a function is a data.frame, the number
#' of rows to display.
#'
#' @export
#'
render_funr <- function(x, max_rows = 100){

  #print(class(out))
  if(class(x)[1] == "try-error")
    cat("")

  out = try(x$value, silent = TRUE)
  vis = ifelse(length(x$visible) == 0, FALSE, x$visible)

  #message("visible status: ", vis)

  if(!vis){
    return(cat(""))

  }else if(is.data.frame(out)){
    message("Showing the first ", max_rows, " rows of the data.frame")
    try(head(out, max_rows), silent = TRUE)

  }else if(class(out) == "help_files_with_topic"){
    ## print help files
    print(out)

  }else if(is.null(out)){
    ## skip NULL
    cat("")

  }else if(is.list(out)){
    ## print list
    print(out)

  }else if(is.atomic(out)){
    cat(out)

  }else if(is.function(out)){
    print(out)

  }else{
    cat("")
  }
}
