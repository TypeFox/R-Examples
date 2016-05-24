print.bayesx <- function(x, model = NULL, ...)
{
  x <- get.model(x, model)
  n <- length(x)
  if(!is.null(model)) {
    start <- stop <- model
    n <- 1L
  } else {
    start <- 1L
    stop <- n
  }
  ncheck <- n > 1L
  nx <- names(x)
  if(is.null(nx))
    nx <- start:stop
  for(i in start:stop) {
    if(ncheck)
      cat("###", nx[i], "\n")
    .print_bayesx(x[[i]])
  }
  if(ncheck) {
    cat("###\n")
    cat("Object consists of", n, "models\n")
  }

  return(invisible(NULL))	
}

