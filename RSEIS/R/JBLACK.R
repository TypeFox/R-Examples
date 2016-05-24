`JBLACK` <-
function(n, acol=rgb(0,0,0))
  {
    if(missing(acol)) {  acol=rgb(0,0,0) } 
    rep(acol, times=n)
  }

