`myGenInv` <-
function(x) 
{
  n <- dim(x)[1]
  nn <-1/n
  return(solve(x+nn)-nn)
}

