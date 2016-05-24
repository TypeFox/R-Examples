.vec <- function(S) {
  Svec = as.matrix(as.vector(S))
  return(Svec)
}

.vech <- function(S) {
  Svech <- as.matrix(S[lower.tri(S, diag = T)])
  return(Svech)
}

.ltri <- function(S) {
  S[upper.tri(S == 1)] = 0
  return(S)
}

.sem_theta_cal <- function(Ldp, Psp, Btp, Php, nup, app, 
                           Ld, Ps, Bt, Ph, nu, ap) {
  theta <- c(.vec(Ld)[.vec(Ldp) != 0], .vech(Ps)[.vech(Psp) != 0],
             .vec(Bt)[.vec(Btp) != 0], .vech(Ph)[.vech(Php) != 0],
             .vec(nu)[.vec(nup) != 0], .vec(ap)[.vec(app) != 0])
  return(theta)
}

.sem_theta_names <- function(Ldp, Psp, Btp, Php, nup, app) {
  if (prod(Ldp == 0) == 1) {
    Ld_names <- NULL
  } else {
    Ldp_idx <- which(Ldp == 1 | Ldp == -1, arr.ind=TRUE)
    Ld_names <- paste("Ld", "[", Ldp_idx[, 1], ",", Ldp_idx[, 2], "]", sep="")    
  }
  if (prod(Psp == 0) == 1) {
    Ps_names <- NULL
  } else {
    Psp_idx <- which(.ltri(Psp) == 1 | .ltri(Psp) == -1, arr.ind=TRUE)
    Ps_names <- paste("Ps", "[", Psp_idx[, 1], ",", Psp_idx[, 2], "]", sep="")    
  }
  if (prod(Btp == 0) == 1) {
    Bt_names <- NULL
  } else {
    Btp_idx <- which(Btp == 1 | Btp == -1, arr.ind=TRUE)
    Bt_names <- paste("Bt", "[", Btp_idx[, 1], ",", Btp_idx[, 2], "]", sep="")    
  }
  if (prod(Php == 0) == 1) {
    Ph_names <- NULL
  } else {
    Php_idx <- which(.ltri(Php) == 1 | .ltri(Php) == -1, arr.ind=TRUE)
    Ph_names <- paste("Ph", "[", Php_idx[, 1], ",", Php_idx[, 2], "]", sep="")    
  }
  if (prod(nup == 0) == 1) {
    nu_names <- NULL
  } else {
    nup_idx <- which(nup == 1 | nup == -1, arr.ind=TRUE)
    nu_names <- paste("nu", "[", nup_idx[, 1], ",", nup_idx[, 2], "]", sep="")  
  }  
  if (prod(app == 0) == 1) {
    ap_names <- NULL
  } else {
    app_idx <- which(app == 1 | app == -1, arr.ind=TRUE)
    ap_names <- paste("ap", "[", app_idx[, 1], ",", app_idx[, 2], "]", sep="")  
  }  
  theta_names <- c(Ld_names, Ps_names, Bt_names, Ph_names, nu_names, ap_names)
  return(theta_names)
}


.sem_implied_moment_cal <- function(Ld, Ps, Bt, Ph, nu, ap) {
  IBtiv <- solve(diag(1, dim(Ld)[2]) - Bt)
  Xi <- IBtiv %*% Ph %*% t(IBtiv)
  Sg <- Ld %*% Xi %*% t(Ld) + Ps
  mu <- nu + Ld %*% IBtiv %*% ap
  implied_moment <- list(Sg = Sg, mu = mu)
  return(implied_moment)
}


.sem_dml_cal <- function(Cyc, ey, Sg, mu) {
  Sgiv <- solve(Sg)
  dml <- - log(det(Sgiv %*% Cyc)) + sum(diag(Sgiv %*% Cyc)) - dim(Sg)[1] + t(ey - mu) %*% Sgiv %*% (ey - mu)
  return(dml)
}


.sem_rpl_cal <- function(theta, thetap, type, gm, dt) {
  if (sum(thetap == -1) > 0) {
    theta_pl <- c(theta[thetap == -1])
    if (type == 'l2') {
      rpl <- gm * sum((theta_pl)^2)
    } else if (type == 'l1') {
      rpl <- gm * sum(abs(theta_pl))
    } else if (type == 'scad') {
      rpl <- sum(gm * (abs(theta_pl) * (abs(theta_pl) <= gm))) + 
        sum((((gm * dt * abs(theta_pl) - 0.5 * (abs(theta_pl)^2 + gm^2))/(dt - 1)) * (abs(theta_pl) > gm & abs(theta_pl) <= (gm * dt)))) + 
        sum((((gm^2) * (dt^2 - 1) / (2 * (dt - 1))) * (abs(theta_pl) > (gm * dt))))
    } else if (type == 'mcp') {
      rpl <- sum(((gm * abs(theta_pl) - (abs(theta_pl)^2) / (2 * dt)) * (abs(theta_pl) <= (gm * dt)))) + 
        sum(((0.5 * dt * gm^2) * (abs(theta_pl) > (gm * dt))))
    } else {rpl <- 0}
  } else {rpl <- 0}
  return(rpl)
}


.sem_threshold <- function(ttq, type, wq, gm, dt) {
  if (type == "l2") {
    ttq <- ttq / (1 + 2 * gm * wq)
  } else if (type == "l1") {
    ttq <- sign(ttq) * max(abs(ttq) - gm * wq, 0)
  } else if (type == "scad") {
    if (abs(ttq) <= gm * (1+wq)) {
      ttq <- sign(ttq) * max(abs(ttq) - gm * wq, 0)
    } else if (abs(ttq) > gm * (1 + wq) & abs(ttq) <= gm * dt) {
      ttq <- sign(ttq) * max(abs(ttq) - gm * wq * dt / (dt - 1), 0) / (1 - (wq / (dt - 1)))
    } else {}
  } else if (type == "mcp") {
    if (abs(ttq) <= gm * dt) {
      ttq <- sign(ttq) * max(abs(ttq) - gm * wq, 0) / (1 - (wq / dt))
    } else {}
  } else {}
  return(ttq)
}

.sem_estep <- function(Cyc, ey, Sg, mu, Ld, Ps, Bt, Ph, nu, ap) {
  Sgiv <- solve(Sg)
  IBtiv <- solve(diag(1, dim(Ld)[2]) - Bt)
  Xi <- IBtiv%*%Ph%*%t(IBtiv)
  KM <- Xi%*%t(Ld)%*%Sgiv
  JM <- IBtiv%*%ap - KM%*%mu
  eet <- JM + KM%*%ey
  Cyet <- ey%*%t(JM) + (Cyc + tcrossprod(ey)) %*% t(KM)
  Cet <- Xi - Xi %*% t(Ld) %*% t(KM) + JM %*% t(JM) + JM %*% t(ey) %*% t(KM) + KM %*% ey %*% t(JM) + 
    KM %*% (Cyc + tcrossprod(ey)) %*% t(KM)
  mis_moment <- list(Cet = Cet, Cyet = Cyet, eet = eet)
  return(mis_moment)
}



.sem_mstep_Ld <- function(Ldp, Ld, Ps, nu, Cyet, Cet, ey, eet, type, gm, dt, P, Psp, Psp_type) {
  if (all(Ld == 0)) {
    Ld <- Ld
  } else {
    if (Psp_type == 1) {
      Psiv <- matrix(0, P, P)
    } else if (Psp_type == 2) {
      Psiv <- diag(1 / diag(Ps))    
    } else if (Psp_type == 3) {
      Psiv <- solve(Ps)
    } else {
      merror_idc = diag(Psp) !=0
      if (all(merror_idc)) {
        Psiv <- solve(Ps)        
      } else {
        Psiv <- matrix(0, P, P)
        Psiv[merror_idc, merror_idc] <- solve(Ps[merror_idc, merror_idc])          
      }
    }
    for (p in 1:P) {
      idx <- which(Ldp[p,] == 1|Ldp[p,] == -1)
      for (j in idx) {
        CRsum <- 0
        if (Psp_type != 1) {
          for (k in (1:P)[(1:P) != p]) {
            CRsum <- CRsum + (Psiv[p, k] / Psiv[p, p]) * (Cyet[k, j] - nu[k, 1] * eet[j, 1] - Ld[k, ] %*% Cet[, j])
          }          
        }
        ldq <- ((Cyet[p, j] - nu[p, 1] * eet[j, 1] - Ld[p, -j] %*% Cet[j, -j]) + CRsum) / Cet[j, j]
        if (Ldp[p, j] == 1) {
          Ld[p, j] <- ldq 
        } else {
          wq <- 1 / (Psiv[p, p] * Cet[j, j])
          Ld[p,j] <- .sem_threshold(ldq, type, wq, gm, dt) 
        }     
      }
    }
  }
  return(Ld)
}


.sem_mstep_Ps <- function(Psp, Ps, Ld, nu, Cyc, Cyet, Cet, ey, eet, type, gm, dt, P, Psp_type) {
  Cep <- Cyc + tcrossprod(ey) - nu %*% t(ey) - Ld %*% t(Cyet) - ey %*% t(nu) + nu %*% t(nu) + 
    Ld %*% eet %*% t(nu) - Cyet %*% t(Ld) + nu %*% t(eet) %*% t(Ld) + Ld %*% Cet %*% t(Ld)
  if (Psp_type == 1) {
    Ps <- Ps
  } else if (Psp_type == 2) {
    Ps <- diag(diag(Cep))    
  } else if (Psp_type == 3) {
    Ps <- Cep
  } else {
    merror_idc <- diag(Psp) !=0
    for (p in 1:P) {
      if (all(merror_idc)) {
        Psiv <- solve(Ps)        
      } else {
        Psiv <- matrix(0, P, P)
        Psiv[merror_idc, merror_idc] <- solve(Ps[merror_idc, merror_idc])          
      }
      idx <- which(Psp[, p] == 1|Psp[, p] == -1)
      for (j in idx[idx > p]) {
        Psivp <- solve(Ps[-p, -p])
        U <- Psivp %*% matrix(Cep[-p, p], (P - 1), 1)
        V <- Psivp %*% Cep[-p, -p] %*% Psivp
        Vsum <- 0
        for (k in idx) {
          if (k < p) {Vk <- V[(j - 1),k]}
          else if (k == p) {Vk <- 0}
          else if (k == j) {Vk <- 0}
          else {Vk <- V[(j - 1), (k - 1)]}
          Vsum <- Vsum + Ps[k,p] * Vk
        }
        psq <- (U[(j - 1)] - Vsum) / V[(j - 1), (j - 1)]
        if (Psp[p, j] == 1) {
          Ps[p, j] <- psq
          Ps[j, p] <- Ps[p, j]
        } else {
          wq <- 1 / (Psiv[p, p] * V[(j - 1), (j - 1)])
          Ps[p, j] <- .sem_threshold(psq, type, wq, gm, dt) 
          Ps[j, p] <- Ps[p, j]
        }
      }
      Psivp <- solve(Ps[-p, -p])
      U <- Psivp %*% matrix(Cep[-p, p], (P - 1), 1)
      V <- Psivp %*% Cep[-p, -p] %*% Psivp
      Vsum <- 0
      Usum <- 0
      for (k in idx) {
        if (k < p) {Uk <- U[k]}
        else if (k == p) {Uk <- 0}
        else {Uk <- U[(k - 1)]}
        Usum <- Usum + Ps[k,p]*Uk
        for (l in idx) {
          if (l < p & k < p) {Vkl <- V[l, k]}
          else if (k == p | l == p) {Vkl <- 0}
          else if (l < p & k > p) {Vkl <- V[l,(k - 1)]}
          else if (l > p & k < p) {Vkl <- V[(l - 1), k]}
          else {Vkl <- V[(l - 1), (k - 1)]}
          Vsum <- Vsum + Ps[k, p] * Ps[l, p] * Vkl
        }
      }
      Ps[p, p] <- Cep[p, p] - 2*Usum + Vsum + matrix(Ps[p, -p], 1, (P - 1)) %*% solve(Ps[-p, -p]) %*% matrix(Ps[-p, p], (P - 1), 1)
    }    
  }
  return(Ps)
}


.sem_mstep_Bt <- function(Btp, Bt, Ph, ap, Cet, eet, type, gm, dt, M, Php_type) {
  if (all(Bt == 0)) {
    Bt <- Bt
  } else {
    if (Php_type == 1) {
      Phiv <- matrix(0, M, M)
    } else if (Php_type == 2) {
      Phiv <- diag(1 / diag(Ph))
    } else if (Php_type == 3) {
      Phiv <- solve(Ph)
    } else {
      Phiv <- solve(Ph)
    }
    for (m in 1:M) {
      idx <- which(Btp[m, ] == 1 | Btp[m, ] == -1)
      for (j in idx) {
        CRsum <- 0
        if (Php_type != 1) {
          for (k in (1:M)[(1:M) != m]) {
            CRsum <- CRsum + (Phiv[m, k] / Phiv[m, m]) * (Cet[k, j] - ap[k, 1] * eet[j, 1] - Bt[k, ] %*% Cet[, j])
          } 
        }
        btq <- ((Cet[m, j] - ap[m, 1] * eet[j, 1] - Bt[m, -j] %*% Cet[j, -j]) + CRsum) / Cet[j, j]
        if (Btp[m, j] == 1) {
          Bt[m, j] <- btq
        } else {
          wq <- 1 / (Phiv[m, m] * Cet[j, j])
          Bt[m, j] <- .sem_threshold(btq, type, wq, gm, dt)
        }
      }
    }    
  }
  return(Bt)
}  


.sem_mstep_Ph <- function(Php, Ph, Bt, ap, Cet, eet, type, gm, dt, M, Php_type) {
  Czt <- Cet - ap %*% t(eet) - Bt %*% Cet - eet %*% t(ap) + ap %*% t(ap) + Bt %*% eet %*% t(ap) - 
    Cet %*% t(Bt) + ap %*% t(eet) %*% t(Bt) + Bt %*% Cet %*% t(Bt)
  if (Php_type == 1) {
    Ph <- Ph
  } else if (Php_type == 2){
    Ph <- diag(diag(Czt))
  } else if (Php_type == 3) {
    Ph <- Czt
  } else {
    for (m in 1:M) {
      Phiv <- solve(Ph)
      idx <- which(Php[, m] == 1 | Php[, m] == -1)
      for (j in idx[idx > m]) {
        Phivm <- solve(Ph[-m, -m])
        U <- Phivm %*% matrix(Czt[-m, m], (M - 1), 1)
        V <- Phivm %*% Czt[-m,-m] %*% Phivm
        Vsum <- 0
        for (k in idx) {
          if (k < m) {Vk <- V[(j-1), k]}
          else if (k == m) {Vk <- 0}
          else if (k == j) {Vk <- 0}
          else {Vk = V[(j - 1), (k - 1)]}
          Vsum <- Vsum + Ph[k, m] * Vk
        }
        phq <- (U[(j - 1)] - Vsum) / V[(j - 1), (j - 1)]
        if (Php[m, j] == 1) {
          Ph[m, j] <- Ph[j, m] <- phq
        } else {
          wq <- 1/(Phiv[m, m] * V[(j - 1), (j - 1)])
          Ph[m,j] <- Ph[j,m] <- .sem_threshold(phq, type, wq, gm, dt) 
        }
      }
      Phivm <- solve(Ph[-m,-m])
      U <- Phivm %*% matrix(Czt[-m, m], (M - 1), 1)
      V <- Phivm %*% Czt[-m, -m] %*% Phivm
      Vsum <- 0
      Usum <- 0
      for (k in idx) {
        if (k < m) {Uk <- U[k]}
        else if (k == m) {Uk <- 0}
        else {Uk <- U[(k - 1)]}
        Usum <- Usum + Ph[k, m] * Uk
        for (l in idx) {
          if (l < m & k < m) {Vkl <- V[l, k]}
          else if (k == m | l == m) {Vkl <- 0}
          else if (l < m & k > m) {Vkl <- V[l, (k - 1)]}
          else if (l > m & k < m) {Vkl <- V[(l - 1), k]}
          else {Vkl <- V[(l - 1), (k - 1)]}
          Vsum <- Vsum + Ph[k, m] * Ph[l, m] * Vkl
        }
      }
      Ph[m, m] <- Czt[m, m] - 2 * Usum + Vsum + matrix(Ph[m, -m], 1, (M - 1)) %*% solve(Ph[-m, -m]) %*% matrix(Ph[-m, m], (M - 1), 1)
    }    
  } 
  return(Ph)
}


.sem_mstep_nu <- function(nup, nu, Ld, Ps, ey, eet, type, gm, dt, P, Psp, Psp_type) {
  if (all(nup == 0)) {
    nu <- nu
  } else {
    if (all(nup == 1)) {
      nu <- ey
    } else {
      if (Psp_type == 1) {
        Psiv <- matrix(0, P, P)
      } else if (Psp_type == 2) {
        Psiv <- diag(1 / diag(Ps))
      } else if (Psp_type == 3) {
        Psiv <- solve(Ps)
      } else {
        merror_idc = diag(Psp) !=0
        if (all(merror_idc)) {
          Psiv <- solve(Ps)        
        } else {
          Psiv <- matrix(0, P, P)
          Psiv[merror_idc, merror_idc] <- solve(Ps[merror_idc, merror_idc])          
        }
        idx <- which(nup[,1] == 1 | nup[, 1] == -1)
        for (p in idx) {
          CRsum <- 0
          for (k in (1:P)[(1:P) != p]) {
            CRsum <- CRsum + (Psiv[p, k] / Psiv[p, p]) * (ey[k, 1] - nu[k,1] - Ld[k,] %*% eet)
          }
          nuq <- ((ey[p,1] - Ld[p,] %*% eet) + CRsum)
          if (nup[p, 1] == 1) {
            nu[p, 1] <- nuq
          } else {
            wq <- 1 / (Psiv[p,p])
            nu[p,1] <- .sem_threshold(nuq, type, wq, gm, dt)
          }
        }
      }
    }
  } 
  return(nu)
}


.sem_mstep_ap = function(app, ap, Bt, Ph, eet, type, gm, dt, M, Php_type) {
  if (all(app == 0)) {
    ap <- ap
  } else {
    if (Php_type == 1) {
      Phiv <- matrix(0, M, M)
    } else if (Php_type == 2) {
      Phiv <- diag(1 / diag(Ph))
    } else if (Php_type == 3) {
      Phiv <- solve(Ph)
    } else {
      Phiv <- solve(Ph)
    }
    idx <- which(app[, 1] == 1 | app[, 1] == -1)
    for (m in idx) {
      CRsum <- 0
      for (k in (1:M)[(1:M) != m]) {
        CRsum <- CRsum + (Phiv[m, k] / Phiv[m, m]) * (eet[k, 1] - ap[k, 1] - matrix(Bt[k, ], 1, M) %*% eet)
      }
      apq <- ((eet[m, 1] - matrix(Bt[m, ], 1, M) %*% eet) + CRsum)
      if (app[m, 1] == 1) {
        ap[m, 1] <- apq
      } else {
        wq <- 1 / (Phiv[m, m])
        ap[m, 1] <- .sem_threshold(apq, type, wq, gm, dt) 
      }
    }    
  } 
  return(ap)
}


.sem_ecm <- function(pattern, value, penalty_fit, control, analysis_info, obs_moment) {
  N <- analysis_info$N
  P <- analysis_info$P
  M <- analysis_info$M
  Qall <- analysis_info$Qall
  
  Ldp <- pattern$Ldp
  Psp <- pattern$Psp
  Btp <- pattern$Btp
  Php <- pattern$Php
  nup <- pattern$nup
  app <- pattern$app
  
  if (all(Psp == 0)) {
    Psp_type = 1
  } else if (all(diag(Psp) == 1) & all(Psp[lower.tri(Psp)] == 0)) {
    Psp_type = 2
  } else if (all(Psp == 1)) {
    Psp_type = 3
  } else {
    Psp_type = 4
  }
  
  if (all(Php == 0)) {
    Php_type = 1
  } else if (all(diag(Php) == 1) & all(Php[lower.tri(Php)] == 0)) {
    Php_type = 2
  } else if (all(Php == 1)) {
    Php_type = 3
  } else {
    Php_type = 4
  }
  
  Ld <- value$Ld
  Ps <- value$Ps
  Bt <- value$Bt
  Ph <- value$Ph
  
  if (all(nup == 1) & all(app == 0) & all(value$ap == 0)) {
    intercept_type = 1
    nu <- obs_moment$mu
    ap <- value$ap
  } else {
    intercept_type = 2
    nu <- value$nu
    ap <- value$ap  
  } 
  
  type <- penalty_fit$type
  gm <- penalty_fit$gm 
  dt <- penalty_fit$dt
  
  itmax <- control$itmax
  eps <- control$eps
  
  Cyc <- obs_moment$Sg
  ey <- obs_moment$mu
  
  implied_moment <- .sem_implied_moment_cal(Ld, Ps, Bt, Ph, nu, ap) 
  Sg <- implied_moment$Sg
  mu <- implied_moment$mu
  dpl <- Inf
  
  for (it in 1:itmax) {
    mis_moment <- .sem_estep(Cyc, ey, Sg, mu, Ld, Ps, Bt, Ph, nu, ap)
    Cet <- mis_moment$Cet
    Cyet <- mis_moment$Cyet
    eet <- mis_moment$eet
    Ld <- .sem_mstep_Ld(Ldp, Ld, Ps, nu, Cyet, Cet, ey, eet, type, gm, dt, P, Psp, Psp_type)
    Ps <- .sem_mstep_Ps(Psp, Ps, Ld, nu, Cyc, Cyet, Cet, ey, eet, type, gm, dt, P, Psp_type)
    Bt <- .sem_mstep_Bt(Btp, Bt, Ph, ap, Cet, eet, type, gm, dt, M, Php_type)
    Ph <- .sem_mstep_Ph(Php, Ph, Bt, ap, Cet, eet, type, gm, dt, M, Php_type)
    if (intercept_type == 2) {
      nu <- .sem_mstep_nu(nup, nu, Ld, Ps, ey, eet, type, gm, dt, P, Psp, Psp_type)
      ap <- .sem_mstep_ap(app, ap, Bt, Ph, eet, type, gm, dt, M, Php_type)
    }
    
    theta <- .sem_theta_cal(Ldp, Psp, Btp, Php, nup, app, Ld, Ps, Bt, Ph, nu, ap)
    thetap <- .sem_theta_cal(Ldp, Psp, Btp, Php, nup, app, Ldp, Psp, Btp, Php, nup, app)
    implied_moment <- .sem_implied_moment_cal(Ld, Ps, Bt, Ph, nu, ap) 
    Sg <- implied_moment$Sg
    mu <- implied_moment$mu
    dml <- .sem_dml_cal(Cyc, ey, Sg, mu)
    rpl <- .sem_rpl_cal(theta, thetap, type, gm, dt)
    dpl_new <- dml + rpl 
    if (abs(dpl - dpl_new) < eps) break
    dpl <- dpl_new
  }
  Q <- sum(theta[thetap == -1] != 0) + sum(thetap == 1)
  df <- P * (P + 3) / 2 - Q
  theta_mat <- list(Ld = Ld, Ps = Ps, Bt = Bt, Ph = Ph, nu = nu, ap = ap)
  rst_ecm <- list(dpl = dpl, dml = dml, Q = Q, df = df, it = it, theta = theta, theta_mat = theta_mat, implied_moment = implied_moment)
  return(rst_ecm)
}


#' A Reference Class for Learning a SEM model via penalized likelihood.

#' @field $data 
#' A N x P data frame contains responses of N observations on P observed variables. In the present version,
#' all variables must be numerical and missing value should be set as \code{NA}. Only complete observations
#' will be used for further analysis

#' @field $pattern 
#' A list of matrices to represent the pattern of a specified SEM model. \code{$pattern} contains six matrices 
#' with element either 0, 1, or -1 to represent fixed, free, or penalized parameter, respectively. 
#' The six matrices are \cr
#' \itemize{
#' \item \code{Ldp}: a P x M pattern matrix for factor loadings. \code{Ldp} must be given and no default 
#' value will be generated through method \code{check()}
#' \item \code{Psp}: a P x P pattern matrix for measurement error covariances. In the present framework, 
#' the diagonal element of \code{Psp} can only be set as 1, i.e., measurement error variance must be 
#' freely estimated. The default value of \code{Psp} is the P x P identity matrix.
#' \item \code{Btp}: a M x M pattern matrix for path coefficients. The default value of \code{Psp} is 
#' the P x P zero matrix.
#' \item \code{Php}: a M x M pattern matrix for latent factors or residuals. In the present framework,
#' the diagonal element of Php can only be set as 1, i.e., residual variance must be freely estimated.
#' The default value of \code{Php} is the M x M identity matrix.
#' \item \code{nup}: a P x 1 pattern matrix for intercepts of observed variables. The default value of 
#' \code{nup} is the P x 1 matrix with all elements being one.
#' \item \code{app}: a M x 1 pattern matrix for intercepts of latent factors. The default value of 
#' \code{app} is the M x 1 matrix with all elements being zero.
#' }

#' @field $value 
#' A list of matrices to specify the starting value and fixed value of a specified SEM model.
#' \code{value} also contains six matrices corresponding to \code{pattern}
#' \itemize{
#' \item \code{Ld}: a P x M matrix for the starting value and fixed value of factor loadings.
#' \code{Ld} must be given and no default value will be generated through method \code{check()}
#' \item \code{Ps}: a P x P matrix for the starting value and fixed value of measurement error covariances.
#' \item \code{Bt}: a M x M matrix for the starting value and fixed value of path coefficients.
#' \item \code{Ph}: a M x M matrix for the starting value and fixed value of latent factor or residual covariances.
#' \item \code{nu}: a P x 1 matrix for the starting value and fixed value of intercepts of observed variables.
#' \item \code{ap}: a M x 1 matrix for the starting value of intercepts of latent factors.
#' }
#' An element in a matrix of \code{value} represents a stating value if the corresponding element in matrix of 
#' \code{pattern} is 1 or -1; otherwise, this element represents a fixed value. 

#' @field $penalty
#' A list of vectors to specify penalization related parameters. \code{penalty} contains three elements: 
#' \itemize{
#' \item \code{type}: a string vector to specify the implemented penalty function. Three penalty functions 
#' can be implemented: \code{"l1"}, \code{"scad"}, and \code{"mcp"}. The default penalty is \code{l1}.
#' \item \code{gm_all}: a numeric vector to specify a candidate set of regularization parameter gamma. The 
#' default value is \code{gm_all = seq(0.01, 0.1, 0.01)}.
#' \item \code{dt_all}: a numeric vector to specify a candidate set of delta. The default value is \itemize{
#' \item \code{dt_all = Inf} if \code{type = "l1"}
#' \item \code{dt_all = c(3, 4)} if \code{type = "scad"}
#' \item \code{dt_all = c(2, 3)} if \code{type = "mcp"}
#' }
#' }

#' @field $control
#' A list of numerical values to specify optimization related parameters. \code{control} contains two elements: 
#' \itemize{
#' \item \code{itmax}: a numeric value to specify the maximal number of iterations. The default value is 500.
#' \item \code{eps}: a numeric value to specify the convergence criterion. The default value is 10^-5.
#' }

#' @field $analysis_info
#' A list of numerical values to represent analysis related numbers. \code{analysis_info} includes four elements:
#' \itemize{
#' \item \code{N}: the number of sample size. 
#' \item \code{P}: the number of observed variables.
#' \item \code{M}: the number of latent factors.
#' \item \code{Qall}: the number of total free and penalized parameters.
#' } 
#' Each element value of \code{analysis_info} will be automatically 
#' assigned when \code{learn()} is executed. 

#' @field $obs_moment
#' A list of sample moment created by excuting \code{learn()}.
#' \itemize{
#' \item \code{Sg}: a P x P sample covariance matrix. 
#' \item \code{mu}: a P x 1 sample mean vector. 
#' } 

#' @field $learn_summary
#' A 3-dimensional array containing the overall analysis result created by executing method \code{learn()}.

#' @field $learn_theta
#' A 3-dimensional array containing the parameter estimates created by executing method \code{learn()}. 

#' @field $fit_summary
#' A matrix containing the overall model information and the values of goodness-of-fit indices obtained by 
#' method \code{fit()}.

#' @field $fit_theta
#' A Qall x 1 matrix containing parameter estimates obtained through method \code{fit()}.

#' @field $fit_value
#' A list of estimated parameter matrices obtained by \code{fit()}.
#' \itemize{
#' \item \code{Ld}: a P x M estimated factor loading matrix.
#' \item \code{Ps}: a P x P estimated measurement error covariance matrix.
#' \item \code{Bt}: a M x M estimated path coefficient matrix.
#' \item \code{Ph}: a M x M estimated residual covariance matrix
#' \item \code{nu}: a P x 1 estimated observed variable intercept
#' \item \code{ap}: a M x 1 estimated latent factor intercept
#' } 

#' @field $fit_moment
#' A list of model implied moments obtained by \code{fit()}. 
#' \itemize{
#' \item \code{Sg}: a P x P estimated model implied covariance.
#' \item \code{mu}: a P x 1 estimated model implied mean.
#' } 

#' @export lslSEM
#' @exportClass lslSEM
#' @import reshape2 ggplot2 methods
#' @examples
#' #Example 1: Factor Analysis Model#
#' #create a P x M population factor loading matrix#
#' Ld0 <- diag(1, 4) %x% matrix(c(.8, .65, .75, .65, .8), 5, 1)
#' Ld0[2, 2] = Ld0[7, 3] = Ld0[12, 4] = Ld0[17, 1] = .5
#' Ld0[9, 1] = Ld0[14, 2] = Ld0[19, 3] = Ld0[4, 4] = .5
#' 
#' #create a M x M population factor covariance matrix#
#' Ph0 <- 0.3 * matrix(1, 4, 1) %*% matrix(1, 1, 4) + diag( .7, 4)
#' 
#' #create a P x P population covariance matrix#
#' Sg0 <- Ld0 %*% Ph0 %*% t(Ld0)
#' diag(Sg0) <- 1
#' 
#' #create a P x P population measurement error covariance matrix#
#' Ps0 <- Sg0 - Ld0 %*% Ph0 %*% t(Ld0)
#' 
#' #create a P x M pattern matrix for factor loadings#
#' Ldp <- 1*(Ld0!=0)
#' Ldp[Ld0 < 0.6] = -1
#' Ldp[1, 1] = Ldp[6, 2] = Ldp[11, 3] = Ldp[16, 4] = 0
#' 
#' #create a M x M pattern matrix for factor covariance#
#' Php <- matrix(1, 4, 4)
#' 
#' #specify field pattern, value, and penalty#
#' pattern <- list(Ldp = Ldp, Php = Php)
#' value <- list(Ld = Ld0, Ps = Ps0, Ph = Ph0)
#' penalty <- list(type = "mcp", gm_all = seq(0.03, .12, .03), dt_all = 1.5)
#' 
#' #generate data with N = 400 and P = 20#
#' Z <- matrix(rnorm(400 * 20, 0, 1), 400, 20)
#' Y <- Z %*% eigen(Sg0)$vectors %*% diag(sqrt(eigen(Sg0)$values)) %*% t(eigen(Sg0)$vectors)
#' Y <- as.data.frame(Y)
#' 
#' #create lslSEM object#
#' rc_sem <- lsl:::lslSEM(data = Y, pattern = pattern, value = value, penalty = penalty)
#' 
#' #check the specification through method check()#
#' rc_sem$check()
#' 
#' #obtain the estimates under each pair of gamma and dt through method learn()#
#' rc_sem$learn()
#' 
#' #obtain the final model based on bic through method fit()#
#' rc_sem$fit(criterion = "bic")
#' 
#' #see overall model information and fit indices of final model#
#' rc_sem$fit_summary
#' 
#' #see estimated Ld of final model#
#' rc_sem$fit_value$Ld
#' 
#' #see the plot for likelihood, aic, and bic under the given gamma and delta values#
#' rc_sem$plot_validation()
#'
#'#################################################################################
#'
#' #Example 2: A general SEM Model#
#' #create a P x M population factor loading matrix#
#' Ld0 <- diag(1, 9) %x% matrix(c(.8,.75,.8), 3, 1)
#' 
#' #create a M x M population path coefficients matrix#
#' Bt0 <- matrix(0, 9, 9)
#' Bt0[2, 1] = Bt0[3, 2] = Bt0[5, 4] = Bt0[6, 5] = Bt0[8, 7] = Bt0[9, 8] =.45
#' Bt0[4, 1] = Bt0[5, 2] = Bt0[6, 3] = Bt0[7, 4] = Bt0[8, 5] = Bt0[9, 6]= .55
#' Bt0iv <- solve(diag(1, 9) - Bt0)
#' 
#' #create a M x M population residual covariance matrix#
#' Ph0 <- diag(0, 9)
#' Ph0[1, 1] <- 1
#' for (m in 2:9) {Ph0[m, m] <- (1 - sum((Bt0iv[m,] ^ 2) * diag(Ph0)))}
#' 
#' #create a P x P population measurement error matrix#
#' Ps0 <- diag(c(.36, 0.4375, .36), 27)
#' #create a P x M population covariance matrix#
#' Sg0 <- Ld0 %*% Bt0iv %*% Ph0 %*% t(Bt0iv) %*% t(Ld0) + Ps0
#' 
#' #create a P x M pattern matrix for factor loadings#
#' Ldp <- (Ld0 != 0)
#' Ldp[1, 1] = Ldp[4, 2] = Ldp[7, 3] = Ldp[10, 4] = Ldp[13, 5] = 0
#' Ldp[16, 6] = Ldp[19, 7] = Ldp[22, 8] = Ldp[25, 9] = 0
#' 
#' #create a P x M pattern matrix for path coefficients#
#' Btp <- matrix(0, 9, 9)
#' Btp[lower.tri(Btp)] <- -1
#' 
#' #specify field pattern, value, and penalty#
#' pattern <- list(Ldp = Ldp, Btp = Btp)
#' value <- list(Ld = Ld0, Bt = Bt0)
#' penalty <- list(type = "mcp", gm_all = seq(0.03, .12, .03), dt_all = 2)
#' 
#' #generate data with N = 400 and P = 27#
#' Z <- matrix(rnorm(400 * 27, 0, 1), 400, 27)
#' Y <- Z %*% eigen(Sg0)$vectors %*% diag(sqrt(eigen(Sg0)$values)) %*% t(eigen(Sg0)$vectors)
#' Y <- as.data.frame(Y)
#' 
#' #create lslSEM object#
#' rc_sem <- lslSEM(data = Y, pattern = pattern, value = value, penalty = penalty)
#' 
#' #check the specification through method check()#
#' rc_sem$check()
#' 
#' #obtain the estimates under each pair of gamma and dt through method learn()#
#' rc_sem$learn()
#' 
#' #obtain the final model based on bic through method fit()#
#' rc_sem$fit(criterion = "bic")
#' 
#' #see overall model information and fit indices of final model#
#' rc_sem$fit_summary
#' 
#' #see estimated Bt of final model#
#' rc_sem$fit_value$Bt
#' 
#' #see the solution path parameters in Bt#
#' rc_sem$plot_path(mat_name = "Bt")

lslSEM <- methods::setRefClass(Class = "lslSEM", 
                               fields = list(
                                 data = "data.frame",
                                 pattern = "list",
                                 value = "list",
                                 penalty = "list",
                                 control = "list",
                                 analysis_info = "list",
                                 obs_moment = "list",
                                 learn_summary = "array",
                                 learn_theta = "array",
                                 fit_summary = "matrix",
                                 fit_theta = "matrix",
                                 fit_value = "list",
                                 fit_moment = "list"),
                               
                               methods = list(    
                                 check = function() {
                                   "Method check() checks the correctness of specifying pattern, value, penalty, and control. If possible,
                                   check() also generates default values for the unspecified field elements."
                                   cat("***check for the field data***\n")
                                   if (length(data) == 0) {
                                     stop("field data is not specified.")
                                   } else {
                                     k <- 0
                                     if (any(!sapply(data, is.numeric))) {
                                       cat("error: some variable in field data is not numeric.\n")
                                       k <- k + 1
                                     } 
                                     if (k == 0) {
                                       cat("field data is OK\n")
                                     } else {
                                       cat(paste(k, "warning(s) or error(s) is generated when checking the field data.\n"))
                                     }
                                   }
                                   cat("***check for the field pattern***\n")
                                   if (length(pattern) == 0) {
                                     stop("field pattern is not specified.")                                     
                                   } else {
                                     if (is.null(pattern$Ldp)) {
                                       stop("pattern$Ldp must be specified.")
                                     } else {
                                       k <- 0
                                       if (any(!is.element(pattern$Ldp, c(-1, 0, 1)))) {
                                         cat("error: element in pattern$Ldp must be -1, 0, or 1.\n")
                                         k <- k + 1
                                       }    
                                       if (is.null(pattern$Psp)) {
                                         pattern$Psp <<- diag(1, dim(pattern$Ldp)[1])
                                         cat("warning: pattern$Psp is not specified and its value is substituted by the default.\n")
                                         k <- k + 1
                                       } else {
                                         if (dim(pattern$Psp)[1] != dim(pattern$Ldp)[1] | dim(pattern$Psp)[1] != dim(pattern$Psp)[2]) {
                                           cat("error: size of pattern$Psp is not correct.\n")
                                           k <- k + 1
                                         } 
                                         if (any(!is.element(pattern$Psp, c(-1, 0, 1)))) {
                                           cat("error: element in pattern$Psp must be -1, 0, or 1.\n")
                                           k <- k + 1
                                         }
                                         if (any(is.element(diag(pattern$Psp), -1))) {
                                           cat("error: element in diag(pattern$Psp) must be 0 or 1.\n")
                                           k <- k + 1
                                         }
                                       }
                                       if (is.null(pattern$Btp)) {
                                         pattern$Btp <<- diag(0, dim(pattern$Ldp)[2])
                                         cat("warning: pattern$Btp is not specified and its value is substituted by the default.\n")
                                         k <- k + 1
                                       } else {
                                         if (dim(pattern$Btp)[1] != dim(pattern$Ldp)[2] | dim(pattern$Btp)[1] != dim(pattern$Btp)[2]) {
                                           cat("error: size of pattern$Btp is not correct.\n")
                                           k <- k + 1
                                         } 
                                         if (any(!is.element(pattern$Btp, c(-1, 0, 1)))) {
                                           cat("error: element in pattern$Btp must be -1, 0, or 1.\n")
                                           k <- k + 1
                                         }
                                         if (any(is.element(diag(pattern$Btp), c(-1, 1)))) {
                                           cat("error: element in diag(pattern$Btp) must be 0.\n")
                                           k <- k + 1
                                         }
                                       }
                                       if (is.null(pattern$Php)) {
                                         pattern$Php <<- diag(1, dim(pattern$Ldp)[2])
                                         cat("warning: pattern$Php is not specified and its value is substituted by the default.\n")
                                         k <- k + 1
                                       } else {
                                         if (dim(pattern$Php)[1] != dim(pattern$Ldp)[2] | dim(pattern$Php)[1] != dim(pattern$Php)[2]) {
                                           cat("error: size of pattern$Php is not correct.\n")
                                           k <- k + 1
                                         } 
                                         if (any(!is.element(pattern$Php, c(-1, 0, 1)))) {
                                           cat("error: element in pattern$Php must be -1, 0, or 1.\n")
                                           k <- k + 1
                                         }
                                         if (any(is.element(diag(pattern$Php), c(-1, 0)))) {
                                           cat("error: element in diag(pattern$Php) must be 1.\n")
                                           k <- k + 1
                                         }
                                       }
                                       if (is.null(pattern$nup)) {
                                         pattern$nup <<- matrix(1, dim(pattern$Ldp)[1], 1)
                                         cat("warning: pattern$nup is not specified and its value is substituted by the default.\n")
                                         k <- k + 1
                                       } else {
                                         if (dim(pattern$nup)[1] != dim(pattern$Ldp)[1] | dim(pattern$nup)[2] != 1) {
                                           cat("error: size of pattern$nup is not correct.\n")
                                           k <- k + 1
                                         } 
                                         if (any(!is.element(pattern$nup, c(-1, 0, 1)))) {
                                           cat("error: element in pattern$nup must be -1, 0, or 1.\n")
                                           k <- k + 1
                                         }
                                       }
                                       if (is.null(pattern$app)) {
                                         pattern$app <<- matrix(0, dim(pattern$Ldp)[2], 1)
                                         cat("warning: pattern$app is not specified and its value is substituted by the default.\n")
                                         k <- k + 1
                                       } else {
                                         if (dim(pattern$app)[1] != dim(pattern$Ldp)[2] | dim(pattern$app)[2] != 1) {
                                           cat("error: size of pattern$app is not correct.\n")
                                           k <- k + 1
                                         } 
                                         if (any(!is.element(pattern$app, c(-1, 0, 1)))) {
                                           cat("error: element in pattern$app must be -1, 0, or 1.\n")
                                           k <- k + 1
                                         }
                                       }
                                       if (k == 0) {
                                         cat("field pattern is OK\n")
                                       } else {
                                         cat(paste(k, "warning(s) or error(s) is generated when checking the field pattern.\n"))
                                       }
                                     }
                                   }
                                   cat("***check for the field value***\n")
                                   if (length(value) == 0) {
                                     stop("field value is not specified.")                                     
                                   } else {
                                     if (is.null(value$Ld)) {
                                       stop("value$Ld must be specified.")
                                     } else {
                                       k <- 0
                                       if (any(dim(value$Ld) != dim(pattern$Ldp))) {
                                         cat("error: size of value$Ld is not correct.\n")
                                         k <- k + 1
                                       }
                                       if (!is.numeric(value$Ld)) {
                                         cat("error: value$Ld is not numeric.\n")
                                         k <- k + 1
                                       }   
                                       if (is.null(value$Ps)) {
                                         value$Ps <<- diag(.1, dim(pattern$Ldp)[1])
                                         cat("warning: value$Ps is not specified and its value is substituted by the default.\n")
                                         k <- k + 1
                                       } else {
                                         if (any(dim(value$Ps) != dim(pattern$Psp))) {
                                           cat("error: size of value$Ps is not correct.\n")
                                           k <- k + 1
                                         }
                                         if (!is.numeric(value$Ps)) {
                                           cat("error: value$Ps is not numeric.\n")
                                           k <- k + 1
                                         }   
                                         if (any(diag(value$Ps) < 0)) {
                                           cat("error: diagonal element of value$Ps can not be negative.\n")
                                           k <- k + 1
                                         }    
                                       }
                                       if (is.null(value$Bt)) {
                                         value$Bt <<- abs(pattern$Btp)* 0.1
                                         cat("warning: value$Bt is not specified and its value is substituted by the default.\n")
                                         k <- k + 1
                                       } else {
                                         if (any(dim(value$Bt) != dim(pattern$Btp))) {
                                           cat("error: size of value$Bt is not correct.\n")
                                           k <- k + 1
                                         }
                                         if (!is.numeric(value$Bt)) {
                                           cat("error: value$Bt is not numeric.\n")
                                           k <- k + 1
                                         }   
                                         if (any(diag(value$Bt) != 0)) {
                                           cat("error: diagonal element of value$Bt must be zero.\n")
                                           k <- k + 1
                                         } 
                                       }
                                       if (is.null(value$Ph)) {
                                         value$Ph <<- diag(.1, dim(pattern$Ldp)[2])
                                         cat("warning: value$Ph is not specified and its value is substituted by the default.\n")
                                         k <- k + 1
                                       } else {
                                         if (any(dim(value$Ph) != dim(pattern$Php))) {
                                           cat("error: size of value$Ph is not correct.\n")
                                           k <- k + 1
                                         }
                                         if (!is.numeric(value$Ph)) {
                                           cat("error: value$Ph is not numeric.\n")
                                           k <- k + 1
                                         }  
                                         if (any(diag(value$Ph) <= 0)) {
                                           cat("error: diagonal element of value$Ph must be positive.\n")
                                           k <- k + 1
                                         }  
                                       }
                                       if (is.null(value$nu)) {
                                         value$nu <<- matrix(colMeans(data, na.rm = T), dim(pattern$Ldp)[1], 1)
                                         cat("warning: value$nu is not specified and its value is substituted by the default.\n")
                                         k <- k + 1
                                       } else {
                                         if (any(dim(value$nu) != dim(pattern$nup))) {
                                           cat("error: size of value$nu is not correct.\n")
                                           k <- k + 1
                                         } 
                                         if (!is.numeric(value$nu)) {
                                           cat("error: value$nu is not numeric.\n")
                                           k <- k + 1
                                         }                                     
                                       }
                                       if (is.null(value$ap)) {
                                         value$ap <<- matrix(0, dim(pattern$Ldp)[2], 1)
                                         cat("warning: value$ap is not specified and its value is substituted by the default.\n")
                                         k <- k + 1
                                       } else {
                                         if (any(dim(value$ap) != dim(pattern$app))) {
                                           cat("error: size of value$ap is not correct.\n")
                                           k <- k + 1
                                         } 
                                         if (!is.numeric(value$ap)) {
                                           cat("error: value$ap is not numeric.\n")
                                           k <- k + 1
                                         }    
                                       }
                                       if (k == 0) {
                                         cat("field value is OK\n")
                                       } else {
                                         cat(paste(k, "warning(s) or error(s) is generated when checking the field value.\n"))
                                       }  
                                     }     
                                   }       
                                   cat("***check for the field penalty***\n")
                                   k <- 0
                                   if (is.null(penalty$type)) {
                                     penalty$type <<- "l1"
                                     cat("warning: penalty$type is not specified and its value is substituted by the default.\n")
                                     k <- k + 1
                                   } else {
                                     if (!is.element(penalty$type, c("l1", "scad", "mcp"))) {
                                       cat("error: penalty$type must be 'l1', 'scad', or 'mcp'.\n")
                                       k <- k + 1
                                     }  
                                   }
                                   if (is.null(penalty$gm_all)) {
                                     penalty$type <<- seq(.01, .1, .01)
                                     cat("warning: penalty$gm_all is not specified and its value is substituted by the default.\n")
                                     k <- k + 1
                                   } else {
                                     if (!is.numeric(penalty$gm_all)) {
                                       cat("error: penalty$gm_all is not numeric.\n")
                                       k <- k + 1
                                     }
                                     if (any(penalty$gm_all < 0)) {
                                       cat("error: some element in penalty$gm_all is negative.\n")
                                       k <- k + 1
                                     }                                     
                                   }
                                   if (is.null(penalty$dt_all)) {
                                     if (penalty$type == "l1") {
                                       penalty$dt_all <<- Inf
                                     } else if (penalty$type == "scad") {
                                       penalty$dt_all <<- c(3, 4)
                                     } else if (penalty$type == "mcp") {
                                       penalty$dt_all <<- c(2, 3)
                                     } else {}
                                     cat("warning: penalty$dt_all is not specified and its value is substituted by the default.\n")
                                     k <- k + 1
                                   } else {
                                     if (!is.numeric(penalty$dt_all)) {
                                       cat("error: penalty$dt_all is not numeric.\n")
                                       k <- k + 1
                                     }
                                     if (any(penalty$dt_all < 0)) {
                                       cat("error: some element in penalty$dt_all is negative.\n")
                                       k <- k + 1
                                     }   
                                   }
                                   if (k == 0) {
                                     cat("field penalty is OK\n")
                                   } else {
                                     cat(paste(k, "warning(s) or error(s) is generated when checking the field penalty.\n"))
                                   }
                                   
                                   cat("***check for the field control***\n")
                                   k <- 0
                                   if (is.null(control$itmax)) {
                                     control$itmax <<- 500
                                     cat("warning: control$itmax is not specified and its value is substituted by the default.\n")
                                     k <- k + 1
                                   } else {
                                     if (length(control$itmax) > 1) {
                                       cat("error: length of control$itmax is larger than 1.\n")
                                       k <- k + 1
                                     }  
                                     if (!is.numeric(control$itmax)) {
                                       cat("error: control$itmax is not numeric.\n")
                                       k <- k + 1
                                     }  
                                   }
                                   
                                   if (is.null(control$eps)) {
                                     control$eps <<- 10^-5
                                     cat("warning: control$eps is not specified and its value is substituted by the default.\n")
                                     k <- k + 1
                                   } else {
                                     if (length(control$eps) > 1) {
                                       cat("error: length of control$eps is larger than 1.\n")
                                       k <- k + 1
                                     }  
                                     if (!is.numeric(control$eps)) {
                                       cat("error: control$eps is not numeric.\n")
                                       k <- k + 1
                                     }  
                                   }
                                   if (k == 0) {
                                     cat("field control is OK\n")
                                   } else {
                                     cat(paste(k, "warning(s) or error(s) is generated when checking the field control.\n"))
                                   }           
                                 },
                                 
                                 learn = function() {
                                   "Method learn() calculates PL estimates under each combination of gm (gamma) and dt (delta) in the field gm_all and dt_all respectively. 
                                   The overall model information will be stored in the field learn_summary and the parameter estimate will be stored in 
                                   the field learn_theta."
                                   if (length(data) != 0) {
                                     analysis_info$N <<- sum(complete.cases(data))
                                     obs_moment$Sg <<- cov(data, use = "pairwise.complete.obs")
                                     obs_moment$mu <<- as.matrix(colMeans(data))                             
                                   } else {
                                     if (length(obs_moment$mu) == 0) {
                                       obs_moment$mu <<- matrix(0, dim(obs_moment$Sg)[1], 1)
                                     }
                                   }
                                   analysis_info$P <<- dim(pattern$Ldp)[1]
                                   analysis_info$M <<- dim(pattern$Ldp)[2]
                                   theta <- .sem_theta_cal(Ldp = pattern$Ldp, Psp = pattern$Psp, 
                                                           Btp = pattern$Btp, Php = pattern$Php,
                                                           nup = pattern$nup, app = pattern$app,
                                                           Ld = value$Ld, Ps = value$Ps,
                                                           Bt = value$Bt, Ph = value$Ph, 
                                                           nu = value$nu, ap = value$ap)
                                   analysis_info$Qall <<- length(theta)
                                   A <- length(penalty$gm_all)
                                   B <- length(penalty$dt_all)
                                   penalty_fit <- list(type = penalty$type)
                                   theta_names <- .sem_theta_names(pattern$Ldp, pattern$Psp, pattern$Btp, pattern$Php, pattern$nup, pattern$app)
                                   gm_name_all <- paste("gm=", as.character(penalty$gm_all), sep="")
                                   dt_name_all <- paste("dt=", as.character(penalty$dt_all), sep="")
                                   learn_summary <<- array(NA, c(7, A, B), dimnames = list(c("dpl", "dml", "Q", "df", "it", "aic", "bic"), gm_name_all, dt_name_all))
                                   learn_theta <<- array(NA, c(analysis_info$Qall, A, B), dimnames = list(theta_names, gm_name_all, dt_name_all))
                                   for (a in A:1) {
                                     penalty_fit$gm <- penalty$gm_all[a]
                                     for (b in B:1) {
                                       penalty_fit$dt <- penalty$dt_all[b]  
                                       rst_ecm <- .sem_ecm(pattern, value, penalty_fit, control, analysis_info, obs_moment)
                                       learn_summary[1:5, a, b] <<- c(rst_ecm$dpl, rst_ecm$dml, rst_ecm$Q, rst_ecm$df, rst_ecm$it)
                                       learn_theta[, a, b] <<- rst_ecm$theta
                                     }    
                                   }
                                   learn_summary[6, , ] <<- learn_summary[2, , ] + (2 / analysis_info$N) * learn_summary[3, , ] 
                                   learn_summary[7, , ] <<- learn_summary[2, , ] + (log(analysis_info$N) / analysis_info$N) * learn_summary[3, , ] 
                                 },
                                 
                                 fit = function(criterion) {
                                   "Method fit() fits a SEM model given the criterion for selecting optimal gm (gamma) and dt (delta). 
                                   Argument criterion can be 'dml' (likelihood), 'aic', or 'bic'. The default value of criterion is 'bic'. 
                                   The final model information and goodness-of-fit will be stored in the field fit_summary. The final 
                                   parameter estimate will be stored in the field fit_theta. The estimated parameter matrix will be stored 
                                   in the field fit_value. The model implied moments will be stored in the field fit_moment."
                                   if (missing(criterion)) {
                                     criterion <- "bic"
                                   }
                                   if (is.element(criterion, c("dml", "aic", "bic"))) {
                                     select_idx <- which(as.matrix(learn_summary[criterion, , ] == min(learn_summary[criterion, , ])), arr.ind = T) 
                                     select_idx <- select_idx[dim(select_idx)[1], , drop = F]
                                     gm <- penalty$gm_all[select_idx[1]]
                                     dt <- penalty$dt_all[select_idx[2]]
                                   } else {
                                     stop("argument criterion must be 'dml', 'aic', or 'bic'")
                                   }
                                   
                                   penalty_fit <- list(type = penalty$type, gm = gm, dt = dt)
                                   theta_names <- .sem_theta_names(pattern$Ldp, pattern$Psp, pattern$Btp, pattern$Php, pattern$nup, pattern$app)
                                   rst_ecm <- .sem_ecm(pattern, value, penalty_fit, control, analysis_info, obs_moment)                           
                                   fit_summary <<- matrix(NA, 18, 1)
                                   colnames(fit_summary) <<- "value"
                                   rownames(fit_summary) <<- c("gamma", "delta", "penalized ml discrepancy", "ml discrepancy", 
                                                               "number of non-zero parameters", "degree of freedom", "number of iterations", 
                                                               "lrt statistic", "srmr", "rmsea", "centrality index", "gamma hat", 
                                                               "cfi", "nnfi", "bl89", "rni", "aic", "bic")
                                   fit_theta <<- as.matrix(rst_ecm$theta)
                                   rownames(fit_theta) <<- theta_names
                                   fit_value$Ld <<- rst_ecm$theta_mat$Ld
                                   fit_value$Ps <<- rst_ecm$theta_mat$Ps
                                   fit_value$Bt <<- rst_ecm$theta_mat$Bt
                                   fit_value$Ph <<- rst_ecm$theta_mat$Ph
                                   fit_value$nu <<- rst_ecm$theta_mat$nu
                                   fit_value$ap <<- rst_ecm$theta_mat$ap
                                   fit_moment$Sg <<- rst_ecm$implied_moment$Sg
                                   fit_moment$mu <<- rst_ecm$implied_moment$mu
                                   
                                   N <- analysis_info$N
                                   P <- analysis_info$P
                                   dml <- rst_ecm$dml
                                   Q <- rst_ecm$Q
                                   lrt <- N * dml
                                   df <- rst_ecm$df
                                   dml_b <- .sem_dml_cal(obs_moment$Sg, obs_moment$mu, diag(diag(obs_moment$Sg)), obs_moment$mu) 
                                   df_b <- P *(P + 3)/2 - 2 * P
                                   lrt_b <- N * dml_b
                                   
                                   fit_summary["gamma", ] <<- gm
                                   fit_summary["delta", ] <<- dt
                                   fit_summary["penalized ml discrepancy", ] <<- rst_ecm$dpl
                                   fit_summary["ml discrepancy", ] <<- rst_ecm$dml
                                   fit_summary["number of non-zero parameters", ] <<- rst_ecm$Q
                                   fit_summary["degree of freedom", ] <<- rst_ecm$df
                                   fit_summary["number of iterations", ] <<- rst_ecm$it
                                   fit_summary["lrt statistic", ] <<- lrt
                                   fit_summary["srmr", ] <<- sqrt(sum(.ltri(((fit_moment$Sg - obs_moment$Sg)^2) / tcrossprod(diag(obs_moment$Sg))))/(P * (P + 1) / 2))
                                   fit_summary["rmsea", ] <<- sqrt(max((lrt - df) / (df * N), 0))
                                   fit_summary["centrality index", ] <<- exp(-0.5 * (lrt - df) / N)
                                   fit_summary["gamma hat", ] <<- P / (P + 2 * (lrt - df) / N)
                                   fit_summary["cfi", ] <<- 1 - max(lrt - df, 0) / max(lrt- df, lrt_b - df_b, 0)
                                   fit_summary["nnfi", ] <<- (lrt_b / df_b - lrt / df)/(lrt_b/ df_b - 1)
                                   fit_summary["bl89", ] <<- (lrt_b - lrt) / (lrt_b- df)
                                   fit_summary["rni", ] <<- ((lrt_b - df_b) - (lrt - df))/(lrt_b - df_b)
                                   fit_summary["aic", ] <<- dml + (2 / N) * Q
                                   fit_summary["bic", ] <<- dml + (log(N) / N) * Q
                                 },  
                                 
                                 plot_validation = function() {
                                   "Method plot_validation() draws a plot for likelihood, aic, and bic under values in gm_all and dt_all"
                                   validation <- reshape2::melt(learn_summary[c("dml", "aic", "bic"), , , drop = F])
                                   validation[, 2] <- as.numeric(substring(validation[, 2], 4))
                                   validation[, 3] <- as.numeric(substring(validation[, 3], 4))
                                   colnames(validation) <- c("criterion", "gm", "dt", "value")
                                   validation_plot <- ggplot2::ggplot(validation, ggplot2::aes(gm, value, colour = criterion)) + ggplot2::geom_line(size = 1) + ggplot2::geom_point(size = 3) + 
                                     ggplot2::facet_grid(.~ dt, labeller = function(variable, value){return(paste("delta = ", as.character(value), sep = ""))}) + 
                                     ggplot2::labs(title = "Plot of Likelihood, aic, and bic", x = "gamma", y = "value")
                                   print(validation_plot) 
                                 },
                                 
                                 plot_path = function(mat_name, row_idx, col_idx) {
                                   "Method plot_path() draws a plot of solution path given mat_name, row_idx, and col_idx.
                                   Argument mat_name must be 'Ld', 'Ps', 'Bt', 'Ph', 'nu', or 'ap'."
                                   theta_names_mat <- matrix(unlist(strsplit(row.names(learn_theta), split="[][,]")), dim(learn_theta)[1], 3, byrow=T)
                                   mat_name_idc <- theta_names_mat[, 1] == mat_name
                                   if (missing(row_idx)) {
                                     row_idc <- rep(T, dim(learn_theta)[1])
                                   } else {
                                     row_idc <- is.element(theta_names_mat[, 2], as.character(row_idx))
                                   }
                                   if (missing(col_idx)) {
                                     col_idc <- rep(T, dim(learn_theta)[1])
                                   } else {
                                     col_idc <- is.element(theta_names_mat[, 3], as.character(col_idx))
                                   }
                                   thetap <- .sem_theta_cal(Ldp = pattern$Ldp, Psp = pattern$Psp, 
                                                            Btp = pattern$Btp, Php = pattern$Php, 
                                                            nup = pattern$nup, app = pattern$app,
                                                            Ld = pattern$Ldp, Ps = pattern$Psp, 
                                                            Bt = pattern$Btp, Ph = pattern$Php, 
                                                            nu = pattern$nup, ap = pattern$app)        
                                   theta_names_penalized <- row.names(learn_theta)[thetap == -1]
                                   path <- reshape2::melt(learn_theta[mat_name_idc & row_idc & col_idc, , , drop = F])
                                   path[, 1] <- as.character(path[,1])
                                   path[, 2] <- as.numeric(substring(path[,2],4))
                                   path[, 3] <- as.numeric(substring(path[,3],4))
                                   penalized <- is.element(path[, 1], theta_names_penalized)
                                   path[, 5] <- penalized
                                   colnames(path) <- c("parameter", "gm", "dt", "value", "penalized")
                                   path_plot <- ggplot2::ggplot(path, ggplot2::aes(gm, value, colour = parameter, linetype = penalized)) + ggplot2::geom_line(size = 1) + 
                                     ggplot2::scale_linetype_manual(values = c("FALSE" = "solid", "TRUE" = "dashed")) +
                                     ggplot2::facet_grid(.~ dt, labeller = function(variable, value) {return(paste("delta = ", as.character(value), sep = ""))}) +
                                     ggplot2::labs(title = paste("Solution Paths of Parameters in ", mat_name, sep = ""), x = "gamma", y = "value")
                                   print(path_plot)
                                 }
                               )
)





