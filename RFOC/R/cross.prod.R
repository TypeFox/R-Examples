`cross.prod` <-
function( B,  A)
  {
    c=list()
    c$x = (B$y*A$z) - (A$y*B$z);
    c$y = (A$x*B$z) - (B$x*A$z);
    c$z = (B$x*A$y) - (A$x*B$y);
    return (c);
  }

