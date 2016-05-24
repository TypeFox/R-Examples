expected <-
  function(x){
    if(!is.matrix(x)){stop("Must be a matrix")}
    rtot <- margin.table(x, 1)
    ctot <- margin.table(x, 2)
    tot <- margin.table(x)
    outer(rtot, ctot, "*")/tot
  }
