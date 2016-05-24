"mono.colors" <-
function (n) {

  if(n<2)
    stop("need to ask for at least 2 colors")
  n=n-1
  col <- rev(gray(0:(n) / (n)))
  col[1] = NA
  return(col)
}

