#compute gradient of QuEST function
get_lambda_J <- function(Q)
{
  rep1p <- rep(1,Q$p)
  
  #Jacobian of support endpoints
  sup_J_fn <- function(Q)
  {
    u_vec <- c(t(Q$u_Fbar))
    Dr <- apply(Q$u_Fbar, c(1,2), function(u) sum(Q$tau^2 / (Q$tau-u)^3 ) )
    temp <- (array(1,dim=dim(Q$u_Fbar)) %o% Q$tau) - (Q$u_Fbar %o% rep1p)
    Q$sup_J <- (Q$u_Fbar %o% Q$tau) / (temp^3 * (Dr %o% rep1p) )
    
    Q$endpoints_J <- 1/Q$n * (Q$u_Fbar^2 %o% rep(1,Q$p)) / temp^2
    if (Q$p == Q$n) Q$endpoints_J[1,1,] <- 0
  }

  sup_J_fn(Q)

  #Jacobian of grid points in u-space
  Q$xi_J_list <- lapply(1:Q$numint, function(i) {
    length_temp <- length(Q$xi_list[[i]])
    temp <- (sin( pi/(2*(length_temp + 1)) * (1:length_temp) ) )^2
    return( tcrossprod((1 - temp), Q$sup_J[i,1,]) + tcrossprod(temp, Q$sup_J[i,2,]) )
  })

  den_J_fn <- function(Q)
  {
    #some objects to be used to in following functions
    TAU_list <- lapply(1:Q$numint, function(i) tcrossprod(rep(1,length(Q$xi_list[[i]])), Q$tau) )
    TAU2_list <- lapply(1:Q$numint, function(i) tcrossprod(rep(1,length(Q$xi_list[[i]])), Q$tau^2) )
    ZXI_list <- lapply(1:Q$numint, function(i) tcrossprod(Q$zxi_list[[i]], rep1p ) )
    XI_list <- lapply(1:Q$numint, function(i) tcrossprod(Re(Q$zxi_list[[i]]), rep1p ) )
    YXI_list <- lapply(1:Q$numint, function(i) tcrossprod(Im(Q$zxi_list[[i]]), rep1p ) )
    YXI2_list <- lapply(1:Q$numint, function(i) tcrossprod(Im(Q$zxi_list[[i]])^2, rep1p ) )
    TAU_XI_list <- lapply(1:Q$numint, function(i) TAU_list[[i]] - XI_list[[i]])
    
    #Jacobian of yxi := imag(zxi)
    Q$yxi_J_list <- lapply(1:Q$numint, function(k) {
      NORM <- TAU_XI_list[[k]]^2 + YXI2_list[[k]]; NORM2 <- NORM*NORM
      Nr <- rowSums(TAU2_list[[k]] * TAU_XI_list[[k]] / NORM2)
      Dr <- rowSums(TAU2_list[[k]] * YXI_list[[k]] / NORM2)
      return ( (TAU_list[[k]] / NORM - TAU2_list[[k]] * TAU_XI_list[[k]] / NORM2 + tcrossprod(Nr, rep1p) * Q$xi_J_list[[k]])  / tcrossprod(Dr, rep1p) )
    })
    
    Q$zxi_J_list <- lapply(1:Q$numint, function(k) Q$xi_J_list[[k]] + 1i*Q$yxi_J_list[[k]])
    
    Q$f_J_list <- lapply(1:Q$numint, function(k) 
      1/(pi*Q$c) * Im(Q$zxi_J_list[[k]] * tcrossprod(1/(Q$zxi_list[[k]]*Q$zxi_list[[k]]), rep1p) ) )
    
    Q$m_LF_J_list <- lapply(1:Q$numint, function(k) {
      TAU_ZXI2 <- (TAU_list[[k]] - ZXI_list[[k]])^2
      1/Q$p * (-ZXI_list[[k]] / TAU_ZXI2 + Q$zxi_J_list[[k]] * tcrossprod(rowSums(TAU_list[[k]] / TAU_ZXI2 ),rep1p) )
    })
    
    Q$x_J_list <- lapply(1:Q$numint, function(k) 
      Re ( Q$zxi_J_list[[k]] * (1 - Q$c*tcrossprod(Q$m_LF_zxi_list[[k]], rep1p) ) - Q$c * ZXI_list[[k]] * Q$m_LF_J_list[[k]] ) )
    
    Q$zeta_J_list <- lapply(1:Q$numint, function(k) 
      1/Q$a * Q$x_J_list[[k]] * tcrossprod(Q$x_list[[k]]^(1/Q$a - 1), rep1p) )
    
    Q$g_J_list <- lapply(1:Q$numint, function(k) 
      Q$a * ( (Q$a - 1) * Q$zeta_J_list[[k]] * tcrossprod(Q$f_list[[k]]*Q$zeta_list[[k]]^(Q$a - 2), rep1p) + 
                tcrossprod(Q$zeta_list[[k]]^(Q$a - 1), rep1p) * Q$f_J_list[[k]] ) )
  }
  
  den_J_fn(Q)
  
  dis_J_fn <- function(Q)
  {
    Q$dis_g_J_list <- lapply(1:Q$numint, function(k) rbind(rep(0,Q$p), Q$g_J_list[[k]], rep(0,Q$p)))
    
    Q$dis_zeta_J_list <- lapply(1:Q$numint, function(k) 
      rbind(1/Q$a * Q$endpoints[k,1]^(1/Q$a - 1) * Q$endpoints_J[k,1,],
            Q$zeta_J_list[[k]],
            1/Q$a * Q$endpoints[k,2]^(1/Q$a - 1) * Q$endpoints_J[k,2,]) )
    if (Q$p == Q$n) Q$dis_zeta_J_list[[1]][1,] <- 0
    
    dis_zeta_J_diff <- lapply(1:Q$numint, function(k) Q$dis_zeta_J_list[[k]][-1,] - Q$dis_zeta_J_list[[k]][-Q$intlength[k],])
    dis_zeta_diff <- lapply(1:Q$numint, function(k) diff(Q$dis_zeta_list[[k]]) )
    dis_g_mean <- lapply(1:Q$numint, function(k) 0.5*(Q$dis_g_list[[k]][-1] + Q$dis_g_list[[k]][-Q$intlength[k]]) )
    dis_g_J_mean <- lapply(1:Q$numint, function(k) 0.5*(Q$dis_g_J_list[[k]][-1,] + Q$dis_g_J_list[[k]][-Q$intlength[k],] ) )
    
    G_J_raw_list <- lapply(1:Q$numint, function(k) rbind(rep(0,Q$p), 
                                              apply(dis_zeta_J_diff[[k]] * tcrossprod(dis_g_mean[[k]], rep1p), 2, cumsum) + 
                                                apply(tcrossprod(dis_zeta_diff[[k]], rep1p) * dis_g_J_mean[[k]], 2, cumsum) ) )
    
    Fends_diff <- Q$Fend - Q$Fstart
    Q$dis_G_J <- lapply(1:Q$numint, function(k)
      Fends_diff[k] / Q$dis_G_raw[[k]][Q$intlength[k]] * (G_J_raw_list[[k]] - 
          tcrossprod(Q$dis_G_raw[[k]], rep1p) * tcrossprod(rep(1,Q$intlength[k]), G_J_raw_list[[k]][Q$intlength[k],]) ) )
  }
  
  dis_J_fn(Q)
  
  lambda_J_fn <- function(Q)
  {
    Q$lambda_J <- matrix(0, nrow=Q$p, ncol = Q$p)
    if ((Q$p - Q$n) <= Q$pzw & Q$pzw > 0) 
        Q$lambda_J[max(1,Q$p - Q$n + 1):Q$pzw, 1:Q$pzw] <- 1 - Q$c
    Q$F_J <- lapply(1:Q$numint, function(k) Q$dis_G_J[[k]][Q$F_idx[[k]],])
    Q$x_J_F <- lapply(1:Q$numint, function(k) 
      tcrossprod((Q$a*Q$dis_zeta_list[[k]][Q$F_idx[[k]]]^(Q$a-1)), rep1p) * Q$dis_zeta_J_list[[k]][Q$F_idx[[k]],] )

    F_J_diff <- lapply(1:Q$numint, function(k) Q$F_J[[k]][-1,] - Q$F_J[[k]][-Q$nidx[k],])
    x_J_F_mean <- lapply(1:Q$numint, function(k) 0.5*(Q$x_J_F[[k]][-1,] + Q$x_J_F[[k]][-Q$nidx[k],]) )
    x_J_F_diff <- lapply(1:Q$numint, function(k) Q$x_J_F[[k]][-1,] - Q$x_J_F[[k]][-Q$nidx[k],])
    X_J_integral <- lapply(1:Q$numint, function(k) 
      F_J_diff[[k]] * tcrossprod(Q$x_F_mean[[k]], rep1p) + tcrossprod(Q$F_diff[[k]], rep1p) * x_J_F_mean[[k]] )
    
    integral_j_kappa <- lapply(1:Q$numint, function(k) {
      q <- Q$quant[[k]]
      b <- Q$bins[[k]]
      tcrossprod(q,rep1p) * Q$x_J_F[[k]][b,] - Q$x_J_F[[k]][b,] * tcrossprod(Q$F[[k]][b], rep1p) -
        Q$F_J[[k]][b,] * tcrossprod(Q$x_F[[k]][b], rep1p) - 
        Q$F_J[[k]][b,] * tcrossprod(((q - Q$F[[k]][b])*Q$x_F_diff[[k]][b] / Q$F_diff[[k]][b]), rep1p) +
        x_J_F_diff[[k]][b,] * tcrossprod( (0.5 * (q - Q$F[[k]][b])^2 / Q$F_diff[[k]][b]), rep1p) -
        F_J_diff[[k]][b,] * tcrossprod( (0.5 * (q - Q$F[[k]][b])^2 * Q$x_F_diff[[k]][b] / Q$F_diff[[k]][b]^2), rep1p) })
    
    X_J_kappa_integral <- lapply(1:Q$numint, function(i) Q$integral_indic[[i]] %*% X_J_integral[[i]] + integral_j_kappa[[i]])
    for (i in 1:Q$numint)
      Q$lambda_J[round(Q$F[[i]][1] * Q$p + 1):round(Q$F[[i]][Q$nidx[i]] * Q$p), ] <- apply(X_J_kappa_integral[[i]], 2, diff) * Q$p
  }
  
  lambda_J_fn(Q)
  return (Q$lambda_J)
}