`in.box` <-
function(x, box, d, boolean=FALSE)
{
  x.box.ind <- rep(TRUE, nrow(x)) 
  for (i in 1:d)
     x.box.ind <- x.box.ind & (box[1,i] <= x[,i]) & (x[,i] <= box[2,i])

  if (boolean)
    return(x.box.ind)
  else  
    return(x[x.box.ind,])
}

