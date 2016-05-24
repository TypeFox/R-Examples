plot_ldlasso <- function( ldlasso.obj ){


  p.vec <- rowSums( abs(ldlasso.obj$cp.obj$beta0.mat) > 1e-6 )
  size.vec <- rowSums( abs(ldlasso.obj$cp.obj$beta0.mat) )

  par(mfrow = c(2,2) )

  plot( ldlasso.obj$bpmap, -ldlasso.obj$log10p, col = 1, pch = 20, type = "p", main = "", xlab = "kB", ylab = "-log10(p)" )
  for( block.bound in ldlasso.obj$block.bounds.vec ){
    lines( rep( block.bound, 2), c(0, max(-ldlasso.obj$log10p) ), lty = 3 )
  }
  
  plot( ldlasso.obj$bpmap, ifelse(abs(ldlasso.obj$beta3)>1e-6,abs(ldlasso.obj$beta3),0), pch = 20, col= 3, type = "p", ylim = c(0,max(abs(c(ldlasso.obj$beta1,ldlasso.obj$beta2,ldlasso.obj$beta3)))), xlab = "kB", ylab = "|beta|" )
  points( ldlasso.obj$bpmap, ifelse(abs(ldlasso.obj$beta2)>1e-6,abs(ldlasso.obj$beta2),0), pch = 20, col = 2 )
  points( ldlasso.obj$bpmap, ifelse(abs(ldlasso.obj$beta1)>1e-6,abs(ldlasso.obj$beta1),0), pch = 20, col = 1 )
  for( block.bound in ldlasso.obj$block.bounds.vec ){
      lines( rep( block.bound, 2), c(0,max(abs(c(ldlasso.obj$beta1,ldlasso.obj$beta2,ldlasso.obj$beta3)))), lty = 3 )
    }

  trace.plot( beta0.mat = ldlasso.obj$cp.obj$beta0.mat, s2.vec = ldlasso.obj$cp.obj$s2.vec, type = "l", indx = NA, s2star = ldlasso.obj$s2star, abs = TRUE )

  plot( ldlasso.obj$cp.obj$s2.vec, ldlasso.obj$cp.obj$cp.vec, type = "l", log = "x", xlab = "s2", ylab = "cp" )
  



}
