`Pradfoc` <-
function(A , MEC, GU, pscale, col)
  {
    C = RPMG::circle()
    
    imageP(MEC$az1, MEC$dip1, MEC$rake1, SCALE=FALSE, UP=MEC$UP, col=col )
    
    lines(C$x, C$y, type='l', col=grey(.6))
    
    
  }

