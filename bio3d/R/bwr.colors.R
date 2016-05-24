"bwr.colors" <-
function (n) {

  # Derived from the function colorpanel
  # by Gregory R. Warnes

  if(n<3)
    warning("not sensible to ask for less than 3 colors")
  
  odd = FALSE
  
  if (n != as.integer(n/2) *2) {
    n   <- n + 1
    odd = TRUE
  }
  low  <- col2rgb("blue")
  mid  <- col2rgb("white")
  high <- col2rgb("red")
  
  lower <- floor(n/2)
  upper <- n - lower
  
  red   <- c(seq(low[1, 1], mid[1, 1], length = lower),
             seq(mid[1,1], high[1, 1], length = upper))/255
  
  green <- c(seq(low[3, 1], mid[3, 1], length = lower),
             seq(mid[3, 1], high[3, 1], length = upper))/255
  
  blue  <- c(seq(low[2, 1], mid[2, 1], length = lower),
             seq(mid[2, 1], high[2, 1], length = upper))/255
    
  if (odd) {
    red   <- red[-(lower + 1)]
    green <- green[-(lower + 1)]
    blue  <- blue[-(lower + 1)]
  }
  rgb(red, blue, green)
}

