draw <- function(x, ...) UseMethod("draw")

draw.vlmc <-
function(x, kind = 3,
         flag = TRUE,
         show.hidden = 0,
         cumulative = TRUE,
         delta = cumulative,
         debug = FALSE, ...)
{
  ## Purpose: Draw a fitted "vlmc" object, see ?vlmc
  ## Author: Martin Maechler, Date: 21 Mar 2000, 12:10

  if(!is.vlmc(x))
      stop("first argument must be a \"vlmc\" object; see ?vlmc")
  ivlmc <- x $ vlmc
  invisible(
  .C(draw_p,
     vlmc.vec     = as.integer(ivlmc),
     size         = length(ivlmc),
     alpha.len    = as.integer  (x$ alpha.len),
     alpha        = as.character(x$ alpha),

     flag         = as.logical(flag),
     debug        = as.logical(debug),
     kind         = as.integer(kind),
     show.hidden  = as.integer(show.hidden),
     cumulative   = as.logical(cumulative),
     delta        = as.logical(delta)))

}
