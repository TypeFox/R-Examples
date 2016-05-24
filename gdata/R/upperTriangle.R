upperTriangle <- function(x, diag=FALSE, byrow=FALSE)
  {
    if(byrow)
      t(x)[rev(upper.tri(x, diag=diag))]
    else
      x[upper.tri(x, diag=diag)]
  }

"upperTriangle<-" <- function(x, diag=FALSE, byrow=FALSE, value)
  {
    if(byrow) {
      ret <- t(x)
      ret[rev(upper.tri(x, diag=diag))] <- value
      t(ret)
    }
    else {        
      x[upper.tri(x, diag=diag)] <- value
      x
    }
  }

lowerTriangle <- function(x, diag=FALSE, byrow=FALSE)
  {
  if(byrow)
    t(x)[rev(lower.tri(x, diag=diag))]
  else
    x[lower.tri(x, diag=diag)]
  }

"lowerTriangle<-" <- function(x, diag=FALSE, byrow=FALSE, value)
  {
  if(byrow) {
    ret <- t(x)
    ret[rev(lower.tri(x, diag=diag))] <- value
    t(ret)
  }
  else {        
    x[lower.tri(x, diag=diag)] <- value
    x
  }
}

