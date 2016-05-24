trace.plot <-
function( beta0.mat, s2.vec, type = "l", indx = NA, s2star = NA, abs = TRUE ){
  if( abs ){
    trace.mat <- abs(t(beta0.mat))
  }else{
    trace.mat <- t(beta0.mat)
  }
  if( abs ){
    plot( 1, type = "n", xlim = range(s2.vec), ylim = range(trace.mat), log = "x", xlab = "s2", ylab = "|beta|" )
  }else{
    plot( 1, type = "n", xlim = range(s2.vec), ylim = range(trace.mat), log = "x", xlab = "s2", ylab = "beta" )
  }
  for( i in 1:nrow(trace.mat) ){
    if( type == "l" ){
      if(!is.na(indx)){
        if( sum( i == indx ) ){
          lines( s2.vec, trace.mat[i,], col = 2 )
        }else{
          lines( s2.vec, trace.mat[i,], col = 1 )
        }
      }else{
        lines( s2.vec, trace.mat[i,], col = 1 )
      }
    }else if ( type == "p" ){
      if(!is.na(indx)){
        if( sum( i == indx ) ){
          points( s2.vec, trace.mat[i,], col = 2 )
        }else{
          points( s2.vec, trace.mat[i,], col = 1 )
        }
      }else{
        points( s2.vec, trace.mat[i,], col = 1 )
      }
    }
  }
  lines( rep(s2star,2), range(trace.mat), lty = 3 )
}

