`dimchecker` <-
function(x,box){

  box.init <- apply(x, 2, range)

  lows <- box[1,]  > box.init[1,]
  
  highs <- box[2,] < box.init[2,]

  rdims <- (lows | highs)

  return(list(either=rdims,lower=lows,upper=highs))

}

