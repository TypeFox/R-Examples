`jpolyval` <-
function(p,x)
  {
    nc = length(p);
    ex = rev(seq(0, nc-1))
    y = rep(0, length(x))
    for( i in 1:length(x))
    y[i] = sum(p * x[i]^ex)
    
    return(y)
  }

