dist.l2d.kgw.u <-
function(x1, varw1, x2, varw2, check=FALSE)  {
  # x1, x2 :       samples
  # varw1, varw2 : bandwidths
  if(check)
  {if(abs(varw1)<.Machine$double.eps | abs(varw2)<.Machine$double.eps)
    {stop("At least one smoothing bandwidth is zero")
    }
  }
  return(sqrt(l2d.kgw.u(x1, varw1, x1, varw1) + l2d.kgw.u(x2, varw2, x2, varw2) - 2*l2d.kgw.u(x1, varw1, x2, varw2)))
}
