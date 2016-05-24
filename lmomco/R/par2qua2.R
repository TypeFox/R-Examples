"par2qua2" <-
function(f,leftpara,rightpara,weight=NULL,...) {
  Q1 <- par2qua(f,leftpara,...)
  Q2 <- par2qua(f,rightpara,...)
  if(is.null(weight)) {
    Q <- (1-f)*Q1 + f*Q2 
  }
  else {
    Q <- (1-weight)*Q1 + weight*Q2
  }
  return(Q)
}

