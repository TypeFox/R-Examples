print.diffusion <-
  function (x, ...) 
  {
    if (!is.diffusion(x)) 
      stop(paste(sQuote("x"), "is not of class", shQuote("diffusion")))
    dum <- deparse(substitute(x))
    cat("\n",shQuote(dum), "is an object of class", shQuote("diffusion"),"which represents")
    cat("\n  a diffusion process, X(t), with:")
    cat("\n   Infinitesimal mean, A1(x,t) =", x$mean)
    cat("\n   Infinitesimal variance, A2(x,t) =", x$var)
    cat("\n\n")
  }