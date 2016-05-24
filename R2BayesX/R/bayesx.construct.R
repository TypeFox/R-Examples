bayesx.construct <- function(object, dir, prg, data) 
{
  UseMethod("bayesx.construct")
}

bayesx.construct.default <- function(object, dir, prg, data) 
{
  cl <- class(object)
  bs <- gsub(".smooth.spec", "", cl, fixed = TRUE)
  info <- if(cl == paste(bs, "smooth.spec", sep = ".")) {
    paste("with basis", sQuote(bs))
  } else {
    paste("of class", sQuote(cl))
  }
  stop(paste("BayesX does not support smooth terms ", info,
    ", it is recommended to use sx() for specifying smooth terms",
    sep = ""))
}
