make.label <-
function(cx, xnam, dimx, vx)
{
  if(cx == "random.bayesx")
    label <- paste("sx(", xnam, sep = "")
  if(cx == "sm.bayesx") {
    if(dimx > 1L)
      label <- paste("sx(", xnam[1L], ",", xnam[2L], sep = "")
    else
      label <- paste("sx(", xnam, sep = "")
  }
  if(cx == "mrf.bayesx")
    label <- paste("sx(", xnam, sep = "")
  if(cx == "geo.bayesx")
    label <- paste("sx(", xnam[1], sep = "")
  if(is.null(vx))
    label <- paste(label, ")", sep = "")
  else
    label <- paste(label, "):", vx, sep = "")

  return(label)
}

