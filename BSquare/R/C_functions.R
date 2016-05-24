################### Matrix Multiplication ########################
# returns A %*% B
# if T_A == 1 => use transpose of A
# if T_B == 1 => use transpose of B
MM_C<-function(A, B, T_A, T_B){
  out <- .C("MM",
            T_A = as.integer(T_A),
            T_B = as.integer(T_B),
            M = as.integer(ifelse(T_A,ncol(A),nrow(A))),            
            N = as.integer(ifelse(T_B,nrow(B),ncol(B))),
            K = as.integer(ifelse(T_A,nrow(A),ncol(A))),
            alpha = as.double(1),
            A = as.double(A),
            B = as.double(B),
            beta = as.double(0),
            C = as.double(rep(0,ifelse(T_A,ncol(A),nrow(A)) * ifelse(T_B,nrow(B),ncol(B))  ))
            )            
  return(out)
}


#make the basis functions
make_B_C<-function(L, kappa, indicator, shape){
  out <- .C("make_B_check",
            L = as.integer(L),
            kappa = as.double(kappa),
            indicator = as.integer(indicator),            
            shape = as.double(shape),
            B = as.double(rep(0,(L+1)*L))
            )            
  return(out)
}

################### Prior Functions #############################
log_det_AR_C<-function(rho,G){
  out <- .C("log_det_AR",
            rho = as.double(rho),
            G=as.integer(G),
            logdet = as.double(0)
            )            
  return(out)
}

sumprec_C<-function(rho,G){
  out <- .C("sum_prec", 
            rho=as.double(rho),
            G=as.integer(G),
            sumprec=as.double(0)
            )            
  return(out)
}

#remember the indexing changes from R to C
sumalpha_C<-function(rho,L,P,p,alpha){
  out <- .C("sum_alpha", 
            rho=as.double(rho),
            L=as.integer(L),
            P=as.integer(P),
            p=as.integer(p),
            alpha=as.double(alpha),
            sumalpha=as.double(0)
            )            
  return(out)
}

#mu is a scalar here, a vector in C code
rhoquad_C<-function(rho, L, P, p, alpha, mu){
  out <- .C("rho_quad",
            rho = as.double(rho),
            L = as.integer(L),
            P = as.integer(P),
            p = as.integer(p), 
            alpha = as.double(alpha), 
            mu = as.double(rep(mu,P)), 
            rhoquad = as.double(0)
            )            
  return(out)
}

################### Quantile Functions #############################

q_norm_C<-function(tau, mn, scale, shape){
  out <- .C("q_norm",
            tau = as.double(tau),
            mn = as.double(mn),
            scale = as.double(scale),
            shape = as.double(shape),
            Q_tau = as.double(0)
            )            
  return(out)
}

q_t_C<-function(tau, mn, scale, shape){
  out <- .C("q_t",
            tau = as.double(tau),
            mn = as.double(mn),
            scale = as.double(scale),
            shape = as.double(shape),
            Q_tau = as.double(0)
            )            
  return(out)
}

q_logistic_C<-function(tau, mn, scale, shape){
  out <- .C("q_logistic",
            tau = as.double(tau),
            mn = as.double(mn),
            scale = as.double(scale),
            shape = as.double(shape),
            Q_tau = as.double(0)
            )            
  return(out)
}

q_alap_C<-function(tau, mn, scale, shape){
  out <- .C("q_alap",
            tau = as.double(tau),
            mn = as.double(mn),
            scale = as.double(scale),
            shape = as.double(shape),
            Q_tau = as.double(0)
            )            
  return(out)
}

q_weibull_C<-function(tau, mn, scale, shape){
  out <- .C("q_weibull",
            tau = as.double(tau),
            mn = as.double(mn),
            scale = as.double(scale),
            shape = as.double(shape),
            Q_tau = as.double(0)
            )            
  return(out)
}

q_gamma_C<-function(tau, mn, scale, shape){
  out <- .C("q_gamma",
            tau = as.double(tau),
            mn = as.double(mn),
            scale = as.double(scale),
            shape = as.double(shape),
            Q_tau = as.double(0)
            )            
  return(out)
}

############# CDFs #####################

p_norm_C<-function(y, mn, scale, shape){
  out <- .C("p_norm",
            y = as.double(y),
            mn = as.double(mn),
            scale = as.double(scale),
            shape = as.double(shape),
            p_y = as.double(0)
            )            
  return(out)
}

p_t_C<-function(y, mn, scale, shape){
  out <- .C("p_t",
            y = as.double(y),
            mn = as.double(mn),
            scale = as.double(scale),
            shape = as.double(shape),
            p_y = as.double(0)
            )            
  return(out)
}

p_logistic_C<-function(y, mn, scale, shape){
  out <- .C("p_logistic",
            y = as.double(y),
            mn = as.double(mn),
            scale = as.double(scale),
            shape = as.double(shape),
            p_y = as.double(0)
            )            
  return(out)
}

p_alap_C<-function(y, mn, scale, shape){
  out <- .C("p_alap",
            y = as.double(y),
            mn = as.double(mn),
            scale = as.double(scale),
            shape = as.double(shape),
            p_y = as.double(0)
            )            
  return(out)
}

p_weibull_C<-function(y, mn, scale, shape){
  out <- .C("p_weibull",
            y = as.double(y),
            mn = as.double(mn),
            scale = as.double(scale),
            shape = as.double(shape),
            p_y = as.double(0)
            )            
  return(out)
}

p_gamma_C<-function(y, mn, scale, shape){
  out <- .C("p_gamma",
            y = as.double(y),
            mn = as.double(mn),
            scale = as.double(scale),
            shape = as.double(shape),
            p_y = as.double(0)
            )            
  return(out)
}


##############log CDFs ###################


log_p_norm_C<-function(y, mn, scale, shape){
  out <- .C("log_p_norm",
            y = as.double(y),
            mn = as.double(mn),
            scale = as.double(scale),
            shape = as.double(shape),
            p_y = as.double(0)
            )            
  return(out)
}

log_p_t_C<-function(y, mn, scale, shape){
  out <- .C("log_p_t",
            y = as.double(y),
            mn = as.double(mn),
            scale = as.double(scale),
            shape = as.double(shape),
            p_y = as.double(0)
            )            
  return(out)
}

log_p_logistic_C<-function(y, mn, scale, shape){
  out <- .C("log_p_logistic",
            y = as.double(y),
            mn = as.double(mn),
            scale = as.double(scale),
            shape = as.double(shape),
            p_y = as.double(0)
            )            
  return(out)
}

log_p_alap_C<-function(y, mn, scale, shape){
  out <- .C("log_p_alap",
            y = as.double(y),
            mn = as.double(mn),
            scale = as.double(scale),
            shape = as.double(shape),
            p_y = as.double(0)
            )            
  return(out)
}
log_p_weibull_C<-function(y, mn, scale, shape){
  out <- .C("log_p_weibull",
            y = as.double(y),
            mn = as.double(mn),
            scale = as.double(scale),
            shape = as.double(shape),
            p_y = as.double(0)
            )            
  return(out)
}

log_p_gamma_C<-function(y, mn, scale, shape){
  out <- .C("log_p_gamma",
            y = as.double(y),
            mn = as.double(mn),
            scale = as.double(scale),
            shape = as.double(shape),
            p_y = as.double(0)
            )            
  return(out)
}

############# densities ################
log_d_norm_C<-function(y, mn, scale, shape){
  out<- .C("log_d_norm",
           y = as.double(y),
           mn = as.double(mn),
           scale = as.double(scale),
           shape = as.double(shape),
           d_Y = as.double(0)
          )
  return(out)
}

log_d_t_C<-function(y, mn, scale, shape){
  out<- .C("log_d_t",
           y = as.double(y),
           mn = as.double(mn),
           scale = as.double(scale),
           shape = as.double(shape),
           d_Y = as.double(0)
           )
  return(out)
}

log_d_alap_C<-function(y, mn, scale, shape){
  out<- .C("log_d_alap",
           y = as.double(y),
           mn = as.double(mn),
           scale = as.double(scale),
           shape = as.double(shape),
           d_Y = as.double(0)
           )
  return(out)
}

log_d_logistic_C<-function(y, mn, scale, shape){
  out<- .C("log_d_logistic",
           y = as.double(y),
           mn = as.double(mn),
           scale = as.double(scale),
           shape = as.double(shape),
           d_Y = as.double(0)
           )
  return(out)
}

log_d_weibull_C<-function(y, mn, scale, shape){
  out<- .C("log_d_weibull",
           y = as.double(y),
           mn = as.double(mn),
           scale = as.double(scale),
           shape = as.double(shape),
           d_Y = as.double(0)
           )
  return(out)
}

log_d_gamma_C<-function(y, mn, scale, shape){
  out<- .C("log_d_gamma",
           y = as.double(y),
           mn = as.double(mn),
           scale = as.double(scale),
           shape = as.double(shape),
           d_Y = as.double(0)
           )
  return(out)
}

threshold_C<-function(L, P, l, alphastar, alpha){
  out<- .C("threshold",
           L = as.integer(L),
           P = as.integer(P),
           l = as.integer(l),
           alphastar = as.double(alphastar),
           alpha = as.double(alpha)
           )
  return(out)
}

clike_C<- function(N, n1, n2, kappa, P, L,
                   X, y, y_low, y_high,
                   beta, alpha,
                   shape, 
                   basis_ind, discrete_ind){
  
  B <- make_B_C(L, kappa, basis_ind, shape)$B
  
    out<- .C("clike_wrapper",
             N = as.integer(N),
             n1 = as.integer(n1),
             n2 = as.integer(n2),
             kappa = as.double(kappa),
             P = as.integer(P),
             L = as.integer(L),
             X = as.double(X),
             y = as.double(y),
             y_low = as.double(y_low),
             y_high = as.double(y_high),
             beta = as.double(beta),
             alpha = as.double(alpha),
             shape = as.double(shape),
             B = as.double(B),
             bin = as.integer(rep(3,N)),
             bin_low = as.integer(rep(3,N)),
             bin_high = as.integer(rep(3,N)),
             ll_sum = as.double(0),
             log_like = as.double(rep(0,N)),
             basis_ind = as.integer(basis_ind),
             discrete_ind = as.integer(discrete_ind)
             )
  return(out)  
  }

MCMC_C<-function(burn, sweeps, 
                 beta, alphastar, shape,
                 X, y, y_low, y_high, 
                 tuning_alpha, tuning_beta, tuning_shape, 
                 beta_eps, alpha_eps,
                 base, L, P, P1, N, n1, n2, n3, n4,
                 mu, sigma2 = rep(1,ncol(X)), rho, 
                 beta_var, shape_var, mu_var, sig_a, sig_b,
                 verbose
                 ){
    basis_ind <- 0;
    if(base == "t"){
      basis_ind <- 1;  
    }
    if(base == "logistic"){
      basis_ind <- 2;
    }
    if(base == "ALAP"){
      basis_ind <- 3;
    }
    if(base == "weibull"){
      basis_ind <- 4;
    }
    if(base == "gamma"){
      basis_ind <- 5;
    }
    out<- .C("MCMC",
             burn = as.integer(burn),
             sweeps = as.integer(sweeps),
             tuning_alpha = as.double(tuning_alpha),
             tuning_beta = as.double(tuning_beta),
             tuning_shape = as.double(tuning_shape),
             beta_eps = as.double(beta_eps),
             alpha_eps = as.double(alpha_eps),
             basis_ind = as.integer(basis_ind), 
             L = as.integer(L),
             N = as.integer(N),
             n1 = as.integer(n1),
             n2 = as.integer(n2),
             n3 = as.integer(n3),
             n4 = as.integer(n4),
             P = as.integer(P),
             P1 = as.integer(P1),
             X = as.double(X),
             y = as.double(y),
             y_low = as.double(y_low),
             y_high = as.double(y_high),
             beta = as.double(beta),
             alphastar = as.double(alphastar),
             shape = as.double(shape),
             mu = as.double(mu),
             sigma2 = as.double(sigma2),
             rho = as.double(rho),
             beta_var = as.double(beta_var),
             shape_var = as.double(shape_var),
             mu_var = as.double(mu_var),
             sig_a = as.double(sig_a),
             sig_b = as.double(sig_b),
             BETA = as.double(rep(0,sweeps*P)),
             ALPHA = as.double(rep(0,sweeps*L*P)),
             MU = as.double(rep(0,sweeps*P)),
             SIGMA2 = as.double(rep(0,sweeps*P)),
             RHO = as.double(rep(0,sweeps*P)),
             SHAPE = as.double(rep(0,sweeps)),
             LPML = as.double(0),
             ACC_BETA = as.integer(rep(0,P)),
             ACC_ALPHA = as.integer(rep(0,L*P)),
             ACC_RHO = as.integer(rep(0,P)),
             ACC_SHAPE = as.integer(0),
             ATT_BETA = as.integer(rep(0,P)),
             ATT_ALPHA = as.integer(rep(0,L*P)),
             verbose = as.integer(verbose)
             )
  return(out)
}

