find.balance.method <- function(method,
                                ...)
{

  if (is.character(method)){

    if (!method == "classical" &
        !method == "stand.diff"){

      stop("Argument 'method' is no valid.")

    }else{
  
      if (method == "classical")
        bal <- balance.classical(...)

      if (method == "stand.diff")
        bal <- balance.stand.diff(...)

    }
  }else{
    stop("Argument 'method' must be a string.")
  }

  return(bal)
  
}
