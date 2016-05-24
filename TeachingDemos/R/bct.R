"bct" <-
function(y,lambda){

  gm <- exp( mean( log(y) ) )

  if(lambda==0) return( log(y)*gm )
  
  yt <- (y^lambda - 1)/( lambda * gm^(lambda-1) )
  return(yt)
}

