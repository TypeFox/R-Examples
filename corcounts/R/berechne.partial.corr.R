`berechne.partial.corr` <-
function(number,Cin,ntemp, upper.triangle.index, upper.triangle.values) {
  index.in <- upper.triangle.index[number,]

  # berechne pc mit vermindertem und vollem Index (entspricht höchstem k)
  pc.vi <- Cin[index.in[1],index.in[2]]
  if (ntemp-sum(index.in==0)>2) {
    for (k in 3:(ntemp-sum(index.in==0))) {

      index.benoetigt.1 <- c(index.in[k],index.in[1])
      if (k>3) {
        for (m in 3:(k-1)) {
          index.benoetigt.1 <- c(index.benoetigt.1, index.in[m])
        }
      }
      index.benoetigt.1 <- c(index.benoetigt.1,0)
      l <- number-1
      while (sum(upper.triangle.index[l,1:k] == index.benoetigt.1)<k) { l <- l-1 }
      pc.1 <- upper.triangle.values[l]

      index.benoetigt.2 <- c(index.in[k],index.in[2])
      if (k>3) {
        for (m in 3:(k-1)) {
          index.benoetigt.2 <- c(index.benoetigt.2, index.in[m])
        }
      }
      index.benoetigt.2 <- c(index.benoetigt.2,0)
      l <- number-1
      while (sum(upper.triangle.index[l,1:k] == index.benoetigt.2)<k) { l <- l-1 }
      pc.2 <- upper.triangle.values[l]

      pc.vi <- (pc.vi-pc.1*pc.2)/(sqrt((1-pc.1^2)*(1-pc.2^2)))
    }
  }

  pc <- pc.vi

  return(pc)
}

