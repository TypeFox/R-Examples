arrow <-
function(x, y, l = .1, w = .3*l, alpha, col.arrow = 'black'){
  w2 <- w/6
  l2 <- l/3
  l3 <- l/2
  x1 <- l*cos(alpha);
  y1 <- l*sin(alpha)
  x2 <- w*cos(alpha+pi/2);
  y2 <- w*sin(alpha+pi/2);
  x3 <- l2*cos(alpha)+w2*cos(alpha+pi/2);
  y3 <- l2*sin(alpha)+w2*sin(alpha+pi/2)
  x4 <- l3*cos(alpha+pi)+w2*cos(alpha+pi/2);
  y4 <- l3*sin(alpha+pi)+w2*sin(alpha+pi/2)
  x5 <- l3*cos(alpha+pi)+w2*cos(alpha+3*pi/2);
  y5 <- l3*sin(alpha+pi)+w2*sin(alpha+3*pi/2)
  x6 <- l2*cos(alpha)+w2*cos(alpha+3*pi/2);
  y6 <- l2*sin(alpha)+w2*sin(alpha+3*pi/2)
  x7 <- w*cos(alpha+3*pi/2);
  y7 <- w*sin(alpha+3*pi/2)

  X <- (par()$usr[2]-par()$usr[1])/par()$pin[1]* c(x1, x2, x3, x4, x5, x6, x7)
  Y <- (par()$usr[4]-par()$usr[3])/par()$pin[2]* c(y1, y2, y3, y4, y5, y6, y7)
  polygon(x + X, y + Y, col = col.arrow, ljoin = 1, border = NA)
}

## Code: Tian, H. and Cazelles, B., \code{WaveletCo}
