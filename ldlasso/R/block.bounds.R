block.bounds <- function( map, block.cood ){
  x.bd <- c()
  wd <- map[length(map)]-map[1]
  for( i in 1:(length(map)+1) ){
    if( i == 1 ){
      x.bd[i] <- map[1]
    }else if( any(i==2:length(map) ) ){
      x.bd[i] <- mean(map[(i-1):i])
    }else if( i==(length(map)+1)){
      x.bd[length(map)+1] <- map[length(map)]
    }
  }
  x.bd[-c(1,length(map)+1)][block.cood[-c(1,length(map)+1)]==1]
}

