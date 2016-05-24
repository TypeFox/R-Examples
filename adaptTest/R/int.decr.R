`int.decr` <-
function (foo, lower, upper, tol=.Machine$double.eps^.5, rel.tol=.Machine$double.eps^.5) {
  if (ge(lower, upper, tol=tol)) {
    result <- 0
    }
  else {
    if (foo(lower) < tol) { 
      result <- 0
      } 
    else { 
      x1 <- lower 
      x2 <- upper 
        if (foo(x2) < tol) { 
          while (abs(x1-x2) >= tol) { 
            x3 <- (x1+x2)/2 
            if (foo(x3) < tol) { 
              x2 <- x3 
              } 
            else { 
              x1 <- x3 
              }             
            } 
          } 
      result <- integrate(foo, lower, x2, rel.tol=rel.tol)$value
      }
    }
  return(result) 
  }

