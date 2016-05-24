message.redef <- function(object, screen.only=F, ...)  { ## by default, two outputs
  ## output to stderr (whether this is screen or not)
  if(is.null(object)) return()
  if(is.list(object)  && length(object)==0) return()
  if(is.data.frame(object)  && nrow(object)==0) return()
  ##ELSE
  messageNAMED(object, ...)
  ## ...stdoutRedirBool: output to stdout  (R_out_....txt file) ()
  if(.blackbox.data$options$stdoutRedirBool &  ! screen.only ) {
    if(is.list(object)) object <- unlist(object)
    if(is.matrix(object)) return(print(object)) ## was: return(invisible(apply(object, 1, cat, fill=T)))
    if(is.numeric(object)) return(print(object)) ## new...
    ## ELSE
    cat(object, "\n")
  }
} ## end redef message
