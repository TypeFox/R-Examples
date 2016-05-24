"getAnisotropy" <- function (model) 
{
  ani <- model@anispec
  if(length(ani$angle) > 0) {
    #non-associative, or associative w/same func. per comp
    if(length(ani$anifunc) == 1) 
      ani$anifunc <- rep(ani$anifunc, 
                         max(model@ncolc))
    model@ncolc <- model@ncolc + length(ani$super)
    ani$calcani <- TRUE
    if(length(ani$useparperp) == 0)
      ani$useparperp <- FALSE		  
  }
  else {
    ani$calcani <- FALSE
    ani$useparperp <- FALSE
  }
  if(length(ani$super) > 0)
    model@cohcol <- model@cohcol + length(ani$super)
  model@anispec <- ani
  model
}


