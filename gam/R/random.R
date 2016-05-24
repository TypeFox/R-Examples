"random" <-
  function(xvar, df = NULL, sigma = 0.)
{
  scall <- deparse(sys.call())
  if(!inherits(xvar, "factor"))
    stop("random() expects a factor or category as its first argument"
         )
  xvar <- C(xvar, rep(0., length(levels(xvar))), 1.)
  attr(xvar, "call") <- substitute(gam.random(data[[scall]], z, w, df = 
                                              df, sigma))
  oldClass(xvar) <- c("smooth", oldClass(xvar))
  xvar
}
