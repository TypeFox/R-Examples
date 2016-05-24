bayesx.construct.generic.smooth.spec <- function(object, dir, prg, data) 
{
  term <- paste(object$term, sep = "", collapse = "*")
  if(!is.null(object$xt)) {
    term <- paste(term, "(", sep = "")
    term <- paste(do.xt(term, object, NULL, TRUE), ")", sep = "")
  }
  if(object$by != "NA")
    term <- make_by(term, object, data)

  return(term)
}

