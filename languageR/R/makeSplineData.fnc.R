`makeSplineData.fnc` <-
function(intr=0) {
  set.seed(314)
  if (!intr) {
    stdev = 2.5
    p = 2
  }
  else {
    stdev = 0.5
    p = 6
  }
  x = seq(0,p*pi,length=100)+2
  if (intr==0) {
    y = 30 + cos(x)
    dfr = data.frame(y=y,X=x)
  } else {
    if (intr==1) {
      z = seq(2,0.2,length=100)
      y = 30 + cos(x)*z
      dfr = data.frame(y=y,X=x,z=z)
    } else {
      z = 1
      y = 30 + cos(x)*z
      dfr = data.frame(y=y,X=x,z=1)
      z = 2
      y = 30 + cos(x)*z
      dfr = rbind(dfr, data.frame(y=y,X=x,z=2))
      z = 3
      y = 30 + cos(x)*z
      dfr = rbind(dfr, data.frame(y=y,X=x,z=3))
      z = 4
      y = 30 + cos(x)*z
      dfr = rbind(dfr, data.frame(y=y,X=x,z=4))
    }
  }
  nsubj = 10
  ranefs = rnorm(10, 0, 3)
  dfr$Subject = rep("S1", nrow(dfr))
  dfr$Ranef = rep(ranefs[1], nrow(dfr)) 
  dfrCopy = dfr
  for (i in 2:10){ 
    tmp = dfrCopy
    tmp$Subject = rep(paste("S", i, sep=""), nrow(tmp))
    tmp$Ranef = ranefs[i]
    dfr = rbind(dfr, tmp)
  }  
  dfr$Subject = as.factor(dfr$Subject)
  dfr$Error = rnorm(nrow(dfr), 0, stdev)
  dfr$Y = dfr$y+dfr$Ranef+dfr$Error
  return(dfr)
}

