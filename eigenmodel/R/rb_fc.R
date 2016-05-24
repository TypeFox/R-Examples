"rb_fc" <-
function( E=Z-ULU(UL) ) {

## sample b from fc
## needs xtx, tx, pm_b and pp_b 

  if(length(pm_b)>0) {
    vr<- solve(  xtx  +  pp_b   )
    mn<- vr%*%( pp_b%*%pm_b +  tx%*%c(E[upper.tri(E)]) )
    b<-rmvnorm(mn,vr)
                      } 
  b 
                   }

