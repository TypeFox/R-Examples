shiftvector <-
function(g,X,x=c(1,0),verbose=FALSE) {
  #
  # compute the optimal shift vector for a calibrated axis.
  #
   costheta <- ((t(g)%*%x))/(sqrt(t(g)%*%g))
   theta <- acos(costheta)
   sintheta <- sin(theta)
   ang <- rad2degree(theta)
   if(verbose) cat("angle vector g -- x-axis:",theta," rad",ang," degrees\n")
   if(g[2]>=0) 
       H <- matrix(c(costheta, sintheta, -sintheta, costheta), ncol = 2) # clockwise
    else
       H <- matrix(c(costheta, -sintheta, sintheta, costheta), ncol = 2) # counterclockwise
   Xr <- X%*%H
   above <- Xr[,2] > 0
   below <- Xr[,2] < 0

   mg <- matrix(rep(1, nrow(X)), ncol = 1)%*%g
   Dp <- diag(as.vector((X%*%g)/(sum(g*g))))
   dmat <- X-Dp%*%mg
   vecdn <- diag(dmat%*%t(dmat) )
   
   if(sum(above)>0) {     
      vecdna <- vecdn[above]
      dmata <- dmat[above,]
      if(is.vector(dmata)) dl <- dmata
      else {
         i <- which.max(vecdna)
         dl <- dmata[i,]
      }
    } else dl <- c(0,0)
   if(sum(below)>0) {
      vecdnb <- vecdn[below]
      dmatb <- dmat[below,]
      if(is.vector(dmatb)) dr <- dmatb
      else {
         i <- which.max(vecdnb)
         dr <- dmatb[i,]
      }
    } else dr <- c(0,0)
   if(verbose) {
     cat("Left shift vector: ",dl,"\n")
     cat("Right shift vector: ",dr,"\n")
   }
   return(list(dr=dr,dl=dl))
}

