seq.ncol <-
function(dt)
  ## auxiliary function to fs.leave.k.out.exact
  {
    numcol <- ncol(dt)
    seq(numcol)
  }
