bayesx.construct.te.smooth.spec <- bayesx.construct.pspline2dimrw2.smooth.spec <-
bayesx.construct.pspline2dimrw1.smooth.spec <- function(object, dir, prg, data)
{
  oby <- if(object$by != "NA") object$by else NULL
  object$by <- object$term[1L]
  object$term <- object$term[2L]
  type <- gsub(".smooth.spec", "", class(object))
  class(object) <- "ps.smooth.spec"
  term <- bayesx.construct(object, dir, prg, data)
  if(type == "pspline2dimrw1") {
    term <- gsub("psplinerw2", "pspline2dimrw1", term)
    term <- gsub("psplinerw1", "pspline2dimrw1", term)
  } else {
    term <- gsub("psplinerw2", "pspline2dimrw2", term)
    term <- gsub("psplinerw1", "pspline2dimrw1", term)
  }
  if(!is.null(oby)) {
    object$by <- oby
    term <- make_by(term, object, data)
  }

  return(term)
}

