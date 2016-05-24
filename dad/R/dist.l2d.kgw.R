dist.l2d.kgw <-
function(x1, varw1, x2, varw2, check = FALSE)  {
  # x1, x2 :       samples
  # varw1, varw2 : bandwidth matrices.
   if(check)
    {if(abs(det(varw1))<.Machine$double.eps | abs(det(varw2))<.Machine$double.eps )
      {stop("One of the smoothing bandwidth matrices is degenerate")
      }
    }  
  return(sqrt(l2d.kgw(x1, varw1, x1, varw1) + l2d.kgw(x2, varw2, x2, varw2) - 2*l2d.kgw(x1, varw1, x2, varw2)))
}
