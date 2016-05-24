## Spearman's rho
rhoFDG <- function(FDGcopula){
  rhoFoo <- dispatch(FDGcopula, "spearman")$rho
  out <- outer(FDGcopula@parameters,FDGcopula@parameters,Vectorize(rhoFoo))
  diag(out) <- rep(1,FDGcopula@dimension)
  out
}
  
## Kendall's tau 
tauFDG <- function(FDGcopula){
  tauFoo <- dispatch(FDGcopula, "kendall")$tau
  out <- outer(FDGcopula@parameters,FDGcopula@parameters,Vectorize(tauFoo))
  diag(out) <- rep(1,FDGcopula@dimension)
  out
}

## Upper tail dependence coefficient 
utdcFDG <- function(FDGcopula){   
  tdcUFoo <- dispatch(FDGcopula, "utdc")$utdc
    out <- outer(FDGcopula@parameters,FDGcopula@parameters,Vectorize(tdcUFoo))
  diag(out) <- rep(1,FDGcopula@dimension)
  out
}
 
## Lower tail dependence coefficient 
ltdcFDG <- function(FDGcopula){   
  tdcLFoo <- dispatch(FDGcopula, "ltdc")$ltdc
  out <- outer(FDGcopula@parameters,FDGcopula@parameters,Vectorize(tdcLFoo))
  diag(out) <- rep(1,FDGcopula@dimension)
  out
}
 
