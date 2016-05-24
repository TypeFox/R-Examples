diffusion <-
  function (text) 
  {
    if (length(text) != 2) 
      stop("object of length 2 is required")
    if (!inherits(text, "list")) 
      text <- as.list(text)
    errors <- !(unlist(lapply(text, is.character)))
    n.errors <- sum(errors)
    if (n.errors > 0) 
      stop(paste("element", rep("s", n.errors > 1), " [", paste((1:2)[errors], 
                                                                collapse = ","), "] in object ", rep(c("is", "are"), 
                                                                                                     c(n.errors == 1, n.errors > 1)), " not character", 
                 sep = ""))
    out <- lapply(text, function(l) as.character(parse(text = l)))
    attr(out, "names") <- c("mean", "var")
    class(out) <- c("diffusion", "list")
    return(out)
  }