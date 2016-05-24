dichromat <- function(colours, type = c("deutan", "protan", "tritan"))
{
  ## transform colors to RGB coordinates
  colours <- t(col2rgb(colours))
  colnames(colours) <- c("r", "g", "b")

  ## compute predicted new dichromat RGB coordinates
  type <- match.arg(type)  
  if(type == "deutan") {
    nred   <- predict(redd,   newdata = colours)
    ngreen <- predict(greend, newdata = colours)
    nblue  <- predict(blued,  newdata = colours)
  } else if(type=="protan") {
    nred   <- predict(redp,   newdata = colours)
    ngreen <- predict(greenp, newdata = colours)
    nblue  <- predict(bluep,  newdata = colours)
  } else if(type=="tritan") {
    nred   <- predict(redt,   newdata = colours)
    ngreen <- predict(greent, newdata = colours)
    nblue  <- predict(bluet,  newdata = colours)
  }

  ## map to unit interval
  nred   <- pmax(0, pmin(1, nred  /255))
  ngreen <- pmax(0, pmin(1, ngreen/255))
  nblue  <- pmax(0, pmin(1, nblue /255))

  ## return color codes
  rgb(nred, ngreen, nblue)
}
