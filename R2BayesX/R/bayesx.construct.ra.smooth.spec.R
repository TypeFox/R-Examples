bayesx.construct.ra.smooth.spec <- bayesx.construct.re.smooth.spec <-
bayesx.construct.random.smooth.spec <- 
function(object, dir, prg, data)
{
  term <- object$term
  if(is.null(object$ins))
    term <- paste(term, "(random", sep = "")
  else
    term <- paste(term, "(hrandom", sep = "")
  term <- paste(do.xt(term, object, NULL), ")", sep = "")
  if(object$by != "NA")
    term <- make_by(term, object, data)

  return(term)
}


bayesx.construct.rsps.smooth.spec <- function(object, dir, prg, data)
{
  term <- paste(object$by, "(hrandom_pspline,centermethod=meansum2", sep = "")
  term <- paste(do.xt(term, object, NULL), ")", sep = "")
  term <- paste(object$term, term , sep = "*")

  return(term)
}
