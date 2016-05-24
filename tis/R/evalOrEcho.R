evalOrEcho <- function(x, resultMode = NULL, n = 0){
  if(is.null(x)) 
    return(NULL)
  if(is.numeric(x)) return(x)
  oldOpt <- options(show.error.messages = FALSE)
  boink <- try(eval.parent(parse(text = x), n = n + 1))
  if(inherits(boink, "try-error") ||
     (!is.null(resultMode) && 
      mode(boink) != resultMode)) 
    boink <- x
  options(oldOpt)
  return(boink)
}
