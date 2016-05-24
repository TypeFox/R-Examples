bayesx.construct.ps.smooth.spec <- bayesx.construct.psplinerw1.smooth.spec <- 
bayesx.construct.psplinerw2.smooth.spec <- bayesx.construct.pspline.smooth.spec <-
function(object, dir, prg, data)
{
  if(length(object$p.order) == 1L) 
    m <- rep(object$p.order, 2L)
  else 
    m <- object$p.order
  m[is.na(m)] <- 2L
  object$p.order <- m
  object$p.order[1L] <- object$p.order[1L] + 1L
  if(class(object) == "psplinerw1.smooth.spec")
    object$p.order[2L] <- 1L
  if(class(object) == "psplinerw2.smooth.spec")
    object$p.order[2L] <- 2L
  if(object$bs.dim < 0L)
    object$bs.dim <- 10L
  if(length(object$p.order) > 1L) {
    if(object$p.order[2L] > 2L) {
      warning("order of the difference penalty not supported by BayesX, set to 2!")
      object$p.order <- c(object$p.order[1L], 2L)
    }
  }
  nrknots <- object$bs.dim - object$p.order[1L] + 1L
  if(nrknots < 5L) {
    warning("number of inner knots smaller than 5 not supported by BayesX, set to 5!",
      call. = FALSE)
    nrknots <- 5L
  }
  termo <- object$term
  term <- paste(termo, "(psplinerw", object$p.order[2L], ",nrknots=",
    nrknots, ",degree=", object$p.order[1L], sep = "")
  object$xt[c("knots", "nrknots", "degree")] <- NULL
  term <- paste(do.xt(term, object, NULL), ")", sep = "")
  if(object$by != "NA")
    term <- make_by(term, object, data)
  
  return(term)
}

