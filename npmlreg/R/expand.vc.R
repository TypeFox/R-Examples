"expand.vc" <-
function(x,ni){
  if (length(ni)==1){
      xx <- x
      for ( i in 2:ni) xx <- rbind(x,xx)
      xx
  } else {
      n <- dim(x)[[1]]
      c <- dim(x)[[2]]
      xx <- NULL
      for ( i in seq(1,n)){
          xx <- rbind(xx,matrix(rep(x[i,],ni[i]),ncol=c,byrow=TRUE))
      }
      xx
  }
}

