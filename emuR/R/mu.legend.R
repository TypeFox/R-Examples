##' make a EMU legend
##' 
##' see function
##' 
##' 
##' @keywords internal
##' @export mu.legend
"mu.legend"<- function(legn, xlim, ylim)
{
  fudge.x <- (xlim[2]-xlim[1])/5
  fudge.y <- (ylim[2]-ylim[1])/5
  if(legn=="tl")
    return(list(x=xlim[1], y=ylim[2]))
  if(legn=="tr")
    return(list(x=xlim[2]-fudge.x, y=ylim[2]))
  if(legn=="br")
    return(list(x=xlim[2]-fudge.x, y=ylim[1]+fudge.y))
  if(legn=="bl")
    return(list(x=xlim[1], y=ylim[1]+fudge.y))
  if(legn=="loc")
    return(graphics::locator(1))
  stop("Unknown legend locator in mu.legend")
}
