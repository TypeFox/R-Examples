summary.Icens <-
function( x, scale=1, ... )
{
  if( attr( x, "model" ) == "MRR" )
    {
      if( is.null( x$rates ) )
        {
          class( x ) <- "glm"
          emat <- ci.lin( x )
        }
      else
        {
          rate.est <- ci.lin( x$rates )
          rate.est[,-(3:4)] <- rate.est[,-(3:4)] * scale
          emat <- rbind( cbind( rate.est, RR=NA )[,c(1:4,7,5:6)],
                        ci.lin( x$cov, Exp=T ) )
        }
    }
  if( attr( x, "model" ) == "AER" )
    {
      rate.est <- ci.lin( x$rates )
      rate.est[,-(3:4)] <- rate.est[,-(3:4)] * scale
      emat <- rate.est
    }
  if( length( x ) == 4 )
    {
      b.est <- x[["boot.ci"]]
      colnames( b.est ) <- c( "Boot-med", "lo", "hi" )
      emat <- cbind( emat, b.est )
    }
  emat
}
