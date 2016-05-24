`INSTFREQS` <-
  function(b,a,w)
  {
    apolyval<-function (coef, z)
      {
        lz <- length(z)
        if (!lz)
          return(numeric(0))
        n <- length(coef)
        if (!n) {
          z[] <- 0
          return(z)
        }
        if (!(mode(coef) == "numeric") && !(mode(coef) == "complex"))
          stop("Argument 'coef' must be a real or complex vector.")
        d_z <- dim(z)
        dim(z) <- lz
        y <- outer(z, (n - 1):0, "^") %*% coef
        dim(y) <- d_z
        return(y)
      }

    
    s = complex(real=0, imaginary=1)*w;
    h = apolyval(b,s) / apolyval(a,s)
    return(h)
  }

