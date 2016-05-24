construct.shrw <-
function(object, dir, prg, data, what)
{
  term <- object$term
  term <- paste(term, "(", what, sep = "")
  term <- paste(do.xt(term, object, NULL), ")", sep = "")
  if(object$by != "NA")
    term <- make_by(term, object, data)

  return(term)
}

