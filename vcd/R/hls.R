hls <- function(h = 1, l = 0.5, s = 1) {
  
  RGB <- function(q1, q2, hue) {
    if (hue > 360) hue <- hue - 360
    if (hue < 0) hue <- hue + 360
    if (hue < 60) 
      q1 + (q2 - q1) * hue / 60
    else if (hue < 180)
      q2
    else if (hue < 240) 
      q1 + (q2 - q1) * (240 - hue) / 60
    else q1
  }

  h <- h * 360
  
  p2 <- if (l <= 0.5)
    l * (1 + s)
  else
    l + s - (l * s)
  p1 <- 2 * l - p2;
  if (s == 0)
    R <- G <- B <- l
  else {
    R <- RGB(p1, p2, h + 120)
    G <- RGB(p1, p2, h)
    B <- RGB(p1, p2, h - 120)
  }
  rgb(R, G, B)
}
