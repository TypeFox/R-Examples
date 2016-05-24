

########################
#III: load C functions # 
########################

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

rhoquad_2_C<-function(rho, G, x1, x2){
  out <- .C("rho_quad_2",
            rho = as.double(rho),
            G = as.integer(G),
            x1 = as.double(x1), 
            x2 = as.double(x2), 
            rhoquad2 = as.double(0)
            )            
  return(out)
}

#returns a sample from a normal (mu, sigma) distribution
rmvnorm_C<- function(mean,sigma){
  out <- .C("rmvnorm", 
            G=as.integer(length(mu)),
            mean=as.double(mean),
            sigma=as.double(sigma),
            y=as.double(rep(0,length(mu)))
            )
  return(out)
}

mspline_knots_C<-function(knots_inter,spline_df){
  knots_inter_length<-length(knots_inter)
  out <- .C("mspline_knots", 
            knots_inter=as.double(knots_inter),
            df_spline=as.integer(spline_df),
            knots_inter_length=as.integer(knots_inter_length),
            TM=as.double(rep(0,2*spline_df + knots_inter_length+ 1))
            )
  return(out)
}

ispline_knots_C<-function(knots_inter,spline_df){
  knots_inter_length<-length(knots_inter)
  out <- .C("ispline_knots", 
            knots_inter=as.double(knots_inter),
            spline_df=as.integer(spline_df),
            knots_inter_length=as.integer(knots_inter_length),
            I_knots=as.double(rep(0,2*spline_df + knots_inter_length+ 2))
            )
  return(out)
}

mspline_C<-function(x,spline_df,i,knots_inter){
  out <- .C("mspline_wrapper", 
            x=as.double(x),
            knots_inter=as.double(knots_inter),
            knots_inter_length=as.integer(length(knots_inter)),
            spline_df=as.integer(spline_df),
            i=as.integer(i-1), #note that I changed the index of i
            v=as.double(0)
            )
  return(out)
}

ispline_C<-function(tau,spline_df,i,knots_inter){
  out <- .C("ispline_wrapper", 
            tau=as.double(tau),
            spline_df=as.integer(spline_df),
            i=as.integer(i-1), #note that I changed the index of i
            knots_inter=as.double(knots_inter),
            knots_inter_length=as.integer(length(knots_inter)),
            v = as.double(0),        
            bin=as.integer(1),
            n = as.integer(1)
            )
  return(out)
}

mspline2_C<-function(tau,spline_df,m,knots_inter){
  out <- .C("mspline2_wrapper", 
            tau=as.double(tau),
            spline_df=as.integer(spline_df),
            m=as.integer(m-1), #note that I changed the index of m
            knots_inter=as.double(knots_inter),
            knots_inter_length=as.integer(length(knots_inter)),
            v = as.double(999),        
            bin=as.integer(1),
            n = as.integer(1)
            )
  return(out)
}

ispline2_C<-function(tau,spline_df,m,knots_inter){
  out <- .C("ispline2_wrapper", 
            tau=as.double(tau),
            spline_df=as.integer(spline_df),
            m=as.integer(m-1), #note that I changed the index of m
            knots_inter=as.double(knots_inter),
            knots_inter_length=as.integer(length(knots_inter)),
            v = as.double(999),        
            bin=as.integer(1),
            n = as.integer(1)
            )
  return(out)
}

mspline3_C<-function(tau,spline_df,m,knots_inter,MKM,M_1){
  out <- .C("mspline3_wrapper", 
            tau=as.double(tau),
            spline_df=as.integer(spline_df),
            m=as.integer(m-1), #note that I changed the index of m
            knots_inter=as.double(knots_inter),
            knots_inter_length=as.integer(length(knots_inter)),
            v = as.double(999),        
            bin=as.integer(1),
            n = as.integer(1),
            MKM = as.double(MKM),
            M_1 = as.integer(M_1)
            )
  return(out)
}

ispline3_C<-function(tau,spline_df,m,knots_inter,IKM,M_1){
  out <- .C("ispline3_wrapper", 
            tau=as.double(tau),
            spline_df=as.integer(spline_df),
            m=as.integer(m-1), #note that I changed the index of m
            knots_inter=as.double(knots_inter),
            knots_inter_length=as.integer(length(knots_inter)),
            v = as.double(999),        
            bin=as.integer(1),
            n = as.integer(1),
            IKM = as.double(IKM),
            M_1 = as.integer(M_1)
            )
  return(out)
}



#returns the leftmost   
findInterval_C<-function(tau,xt){
  out<-.C("find_int",
          n = as.integer(length(xt)), 
          xt = as.double(xt), 
          tau = as.double(tau), 
          ilo = as.integer(0)
          )
  return(out)
}

#function that generates a MVN RV with mean mean mu and cov sigma
rmvnorm_C<-function(mu,Sigma){
  D<-length(mu)
  out <- .C("rmvnorm", 
            y=as.double(rep(0,D)),
            D=as.integer(D),
            mean=as.double(mu),
            sigma=as.double(matrix(diag(D),nrow=D,ncol=D))
            )
  return(out)
}

#function that multiplies 2 matrices
MM_C<-function(M1,M2){
  out <- .C("MM", 
            nrow_M1=as.integer(nrow(M1)),
            ncol_M1=as.integer(ncol(M1)),
            M1=as.double(M1),
            ncol_M2=as.integer(ncol(M2)),
            M2=as.double(M2),
            M3=as.double(matrix(0,nrow=nrow(M1),ncol=ncol(M2)))
            )
  return(out)
}

rootfind_C<-function(M_knots_length, I_knots_length, M, spline_df, tau_scalar, w,
                     y_scalar,M_knots, I_knots, bin, q_low, q_high){
  out <- .C("rootfind", 
            M_knots_length=as.integer(M_knots_length),
            I_knots_length=as.integer(I_knots_length),
            M=as.integer(M),
            spline_df=as.integer(spline_df),
            tau_scalar=as.double(tau_scalar),
            w=as.double(w),
            y_scalar=as.double(y_scalar),
            M_knots=as.double(M_knots),
            I_knots=as.double(I_knots),
            bin=as.integer(bin),
            q_low = as.double(q_low),     
            q_high = as.double(q_high),
            iter_flag = as.integer(0)
            )
  return(out)
}

rootfind_GPU_C<-function(M_knots_length, I_knots_length, M, spline_df, tau_scalar, w,
                     y_scalar,M_knots, I_knots, bin, q_low, q_high, IKM, MKM, M_1, reset_value){
  out <- .C("rootfind_GPU", 
            M_knots_length=as.integer(M_knots_length),
            I_knots_length=as.integer(I_knots_length),
            M=as.integer(M),
            spline_df=as.integer(spline_df),
            tau_scalar=as.double(tau_scalar),
            w=as.double(w),
            y_scalar=as.double(y_scalar),
            M_knots=as.double(M_knots),
            I_knots=as.double(I_knots),
            bin=as.integer(bin),
            q_low = as.double(q_low),     
            q_high = as.double(q_high),
            IKM = as.double(IKM),
            MKM = as.double(MKM),
            M_1 = as.integer(M_1),
            reset_value = as.double(reset_value)
            )
  return(out)
}

q_C<-function(M, tau, w, spline_df, M_knots){
  out <- .C("q_wrapper", 
            M=as.integer(M),
            tau= as.double(tau),     
            w=as.double(w),
            spline_df=as.integer(spline_df),
            M_knots=as.double(M_knots),  
            r_Q = as.double(0)
            )
  return(out)
}

Q_C<-function(M, tau, w, spline_df, I_knots, bin){
  r_Q<-w[1]
  out <- .C("Q_wrapper", 
            M=as.integer(M),
            tau= as.double(tau),     
            w=as.double(w),
            spline_df=as.integer(spline_df),
            I_knots=as.double(I_knots),
            bin=as.integer(bin),     
            r_Q = as.double(r_Q)
            )
  return(out)
}



#build the covariance matrix
#sigma is the standard deviation
#rho is the correlation - must be less than sigma


make_prec_C<-function(rho,G){
  out <- .C("make_prec", 
            rho=as.double(rho),
            G=as.integer(G),
            OMEGA=as.double(rep(0,G^2))
            )            
  return(out)
}

make_prec_2_C<-function(rho,G){
  out <- .C("make_prec_2", 
            rho=as.double(rho),
            G=as.integer(G),
            OMEGA=as.double(rep(0,G^2))
            )            
  return(out)
}

#G is number of rows of square matrix OMEGA

log_det_C<-function(G,OMEGA){
  out <- .C("log_det", 
            G=as.integer(G),
            OMEGA=as.double(OMEGA),
            logdet = as.double(0)
            )            
  return(out)
}

log_det_2_C<-function(rho,G){
  out <- .C("log_det_2",
            rho = as.double(rho),
            G=as.integer(G),
            logdet = as.double(0)
            )            
  return(out)
}

sumtheta_C<-function(rho,M,G,P,m,p,theta){
  out <- .C("sum_theta", 
            rho=as.double(rho),
            M=as.integer(M),
            G=as.integer(G),
            P=as.integer(P),
            m=as.integer(m),
            p=as.integer(p),
            theta=as.double(theta),
            sumtheta=as.double(0)
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

AR_chol_C<-function(rho,G){
  out <- .C("AR_chol", 
            rho=as.double(rho),
            G=as.integer(G),
            G_chol=as.double(rep(0,G^2))
            )            
  return(out)
}

rhoquad_C<-function(rho, M, G, P, m, p, theta, mu){
  out <- .C("rho_quad",
            rho = as.double(rho),
            M = as.integer(M),
            G = as.integer(G),
            P = as.integer(P),
            m = as.integer(m), 
            p = as.integer(p), 
            theta = as.double(theta), 
            mu = as.double(mu), 
            rhoquad = as.double(0)
            )            
  return(out)
}

#function that passes pointer to function 
#(allows me to call different covariances, priors, etc. without renaming in C)
#dyn.load("chris_is_brilliant.dll")
main_function_C<-function(a,b,indicator){
  out <- .C("main_function",
            a= as.double(a),
            b= as.double(b),
            answer = as.double(0),
            indicator=as.integer(indicator)
            )            
  return(out)
}

log_det_AR_C<-function(rho,G){
  out <- .C("log_det_AR",
            rho = as.double(rho),
            G=as.integer(G),
            logdet = as.double(0)
            )            
  return(out)
}

#chol2inv_C returns the lower triangular part of the inverse of the matrix
chol2inv_C<-function(G,OMEGA){
  out <- .C("chol2inv", 
            G=as.integer(G),
            OMEGA=as.double(OMEGA),
            OMEGA1=as.double(rep(0,length(OMEGA)))
            )            
  return(out)
}

#prior function for updating theta 
#PREC is G x G x M x P dimensional array

lp_vague_C<-function(theta){
  out <- .C("lp_vague", 
            theta_length=as.integer(length(theta)),
            theta = as.double(theta),
            lp_sum = as.double(0)
            )            
  return(out)
}

proper_theta_C<-function(M, G, P, m, g, p, theta){
  out <- .C("proper_theta", 
            M=as.integer(M),
            P=as.integer(P),
            G=as.integer(G),
            m=as.integer(m),
            g=as.integer(g),
            p=as.integer(p),
            theta = as.double(theta),
            proper_flag = as.integer(0)
            )            
  return(out)
}  

threshold_C<-function(M, G, P, m, g, thetastar, theta){
  out <- .C("threshold", 
            M=as.integer(M),
            G=as.integer(G),
            P=as.integer(P),
            m=as.integer(m),
            g=as.integer(g),
            thetastar = as.double(thetastar),
            theta = as.double(theta),
            slope_update = as.integer(0)
            )            
  return(out)
}

dotproduct_C<-function(v1,v2){
  out <- .C("dot_product", 
            D_vector=as.integer(length(v1)),
            v1=as.double(v1),
            v2=as.double(v2),
            dp=as.double(0)
            )            
  return(out)
}

#theta length is (P+1)*G*M
#theta is G X (P+1) X M array 
#alpha is G X 1 X M array
#beta is G X P X M array
#rho is 1 X (P+1)
#sigma2 is 1 X (P+1)
#xi_low is 1 X G 
#xi_high is 1 X G
#y is N_g x G matrix
#X is N_g x P x G array
#N 1 X 1 is the number of observations at each gestational age
#tuning_parms is G x P array

proper_theta_C<-function(M, G, P, theta){
  out <- .C("proper_theta", 
            M=as.integer(M),
            G=as.integer(G),
            P=as.integer(P),
            theta = as.double(theta),
            proper_flag = as.integer(0)
            )            
  return(out)
}

gamma_sample_C<-function(shape,scale,G, N){
  
  out <- .C("gamma_sample", 
            N=as.integer(N),
            G=as.integer(G),
            shape=as.double(shape),
            scale = as.double(scale),
            sample= as.double(rep(0,N))
            )
  return(out)
}     

geom_sample_C<-function(N,eps){
    out <- .C("geom_sample", 
            N=as.integer(N),
            eps=as.double(eps),
            sample=as.integer(rep(0,N))
            )
  return(out)
}

#################### MCMC functions ###############################

MCMC_IP_C<-function (burn,iters,
                     tuning_parms, tuning_tail,
                     cbf_eps, theta_eps, 
                     M,  P,  P1, N,  n1, n2, n3, n4, n5,
                     y,  y_low, y_high, X,
                     M_knots_length, I_knots_length, spline_df, 
                     M_knots, I_knots,
                     M_low,  M_high, I_low,  I_high,   
                     thetastar,  mu,  sigma2,  rho,  xi_low,  xi_high,
                     q_low,  q_high, xi_zero,
                     mu_var, cbf_var, tail_mean, tail_var,
                     sig_a, sig_b,
                     IKM, MKM, M_1,verbose){
      out <- .C("MCMC_IP", 
              burn = as.integer(burn),
              iters = as.integer(iters),
              tuning_parms = as.double(tuning_parms),
              tuning_tail = as.double(tuning_tail),
              cbf_eps = as.double(cbf_eps),
              theta_eps = as.double(theta_eps),
              M = as.integer(M),
              P = as.integer(P),
              P1 = as.integer(P1),
              N = as.integer(N),
              n1 = as.integer(n1),
              n2 = as.integer(n2),
              n3 = as.integer(n3),
              n4 = as.integer(n4),
              n5 = as.integer(n5),
              y = as.double(y),
              y_low = as.double(y_low),
              y_high = as.double(y_high),
              X = as.double(X),
              M_knots_length = as.integer(M_knots_length),
              I_knots_length = as.integer(I_knots_length),
              spline_df = as.integer(spline_df),
              M_knots = as.double(M_knots),
              I_knots = as.double(I_knots),
              M_low = as.double(M_low),
              M_high = as.double(M_high),
              I_low = as.double(I_low),
              I_high = as.double(I_high),
              thetastar = as.double(thetastar),
              mu = as.double(mu),
              sigma2 = as.double(sigma2),
              rho = as.double(rho),
              xi_low = as.double(xi_low),
              xi_high = as.double(xi_high),
              q_low = as.double(q_low),
              q_high = as.double(q_high),
              xi_zero = as.integer(xi_zero),
              mu_var = as.double(mu_var),
              cbf_var = as.double(cbf_var),
              tail_mean = as.double(tail_mean),
              tail_var = as.double(tail_var),
              sig_a = as.double(sig_a),
              sig_b = as.double(sig_b),
              THETA = as.double(rep(0,iters*M*P)),
              MU = as.double(rep(0,iters*P)),
              SIGMA2 = as.double(rep(0,iters*P)),
              RHO = as.double(rep(0,iters*P)),
              XI_LOW = as.double(rep(0,iters)),
              XI_HIGH = as.double(rep(0,iters)),
              LPML = as.double(0),
              ACC_THETA = as.integer(rep(0,M*P)),
              ATT_THETA = as.integer(rep(0,M*P)),
              ACC_RHO = as.integer(rep(0,P)),
              ACC_TAIL = as.integer(rep(0,2)),
              IKM = as.double(IKM),
              MKM = as.double(MKM),
              M_1 = as.integer(M_1),
              verbose = as.integer(verbose)
              )
    return(out)
}

MCMC_CP_C<-function (burn,sweeps,
                     tuning_parms, tuning_tail,
                     cbf_eps, theta_eps, 
                     M,  G,  P,
                     M0g, Pg, P1,
                     N,  n1, n2, n3, n4, n5,
                     X, y,  y_low, y_high, 
                     M_knots_length, I_knots_length, spline_df, 
                     M_knots, I_knots,
                     M_low,  M_high, I_low,  I_high,  
                     thetastar,  mu,  sigma2,  rho,  xi_low,  xi_high,
                     q_low,  q_high,
                     xi_zero,
                     mu_var, cbf_var, tail_mean, tail_var,
                     sig_a, sig_b,
                     IKM, MKM, M_1){
      out <- .C("MCMC_CP", 
              burn = as.integer(burn),
              sweeps = as.integer(sweeps),
              tuning_parms = as.double(tuning_parms),
              tuning_tail = as.double(tuning_tail),
              cbf_eps = as.double(cbf_eps),
              theta_eps = as.double(theta_eps),
              M = as.integer(M),
              P = as.integer(P),
              G = as.integer(G),
              N = as.integer(N),
              n1 = as.integer(n1),
              n2 = as.integer(n2),
              X = as.double(X),
              y = as.double(y),
              y_low = as.double(y_low),
              y_high = as.double(y_high),
              M_knots_length = as.integer(M_knots_length),
              I_knots_length = as.integer(I_knots_length),
              spline_df = as.integer(spline_df),
              M_knots = as.double(M_knots),
              I_knots = as.double(I_knots),
              M_low = as.double(M_low),
              M_high = as.double(M_high),
              I_low = as.double(I_low),
              I_high = as.double(I_high),
              thetastar = as.double(thetastar),
              mu = as.double(mu),
              sigma2 = as.double(sigma2),
              rho = as.double(rho),
              xi_low = as.double(xi_low),
              xi_high = as.double(xi_high),
              q_low = as.double(q_low),
              q_high = as.double(q_high),
              xi_zero = as.integer(xi_zero),
              mu_var = as.double(mu_var),
              cbf_var = as.double(cbf_var),
              tail_mean = as.double(tail_mean),
              tail_var = as.double(tail_var),
              sig_a = as.double(sig_a),
              sig_b = as.double(sig_b),
              THETA = as.double(rep(0,sweeps*M*G*P)),
              MU = as.double(rep(0,sweeps*M*P)),
              SIGMA2 = as.double(rep(0,sweeps*M*P)),
              RHO = as.double(rep(0,sweeps*M*P)),
              XI_LOW = as.double(rep(0,sweeps)),
              XI_HIGH = as.double(rep(0,sweeps)),
              TAU = as.double(rep(0,sweeps*N)),
              LPML = as.double(0),
              ACC_THETA = as.integer(rep(0,M*G*P)),
              ATT_THETA = as.integer(rep(0,M*G*P)),
              ACC_RHO = as.integer(rep(0,M*P)),
              ACC_TAIL = as.integer(rep(0,2)),
              IKM = as.double(IKM),
              MKM = as.double(MKM),
              M_1 = as.integer(M_1)    
              )    
      return(out)
}

rhoquad_vector_C<-function(rho, M, G, P, m, p, theta, mu){
  out <- .C("rho_quad_vector",
            rho = as.double(rho),
            M = as.integer(M),
            G = as.integer(G),
            P = as.integer(P),
            m = as.integer(m),
            p = as.integer(p),
            theta = as.double(theta),
            mu = as.double(mu),
            rhoquad2 = as.double(0)
            )            
  return(out)
}


make_betahat_AR_C<-function(G, GG, D, M, m, p, theta, rho){
  out <- .C("make_betahat_AR", 
            G=as.integer(G),
            GG=as.double(GG),
            D=as.integer(D),
            M = as.integer(M),
            m = as.integer(m),
            p = as.integer(p),
            theta = as.double(theta),
            rho = as.double(rho),
            V = as.double(rep(7,D^2)),
            chol_V = as.double(rep(4,D^2)),
            chol_V_inv = as.double(rep(3,D^2)),
            V_inv = as.double(rep(5,D^2)),
            betahat = as.double(rep(5,D))
            )
  return(out)
}

sum_theta_SP_C<-function(rho, M, G, P, g, p, theta){
  out <- .C("sum_theta_SP", 
            rho=as.double(rho),
            M=as.integer(M),
            G = as.integer(G),
            P = as.integer(P),
            g = as.integer(g),
            p = as.integer(p),
            theta = as.double(theta),
            sumtheta = as.double(0)
            )
  return(out)
}

#mu is G x P
#rho is scalar
rho_quad_SP_C<-function(rho, M, G, P, g, p, theta, mu){
  out <- .C("rho_quad_SP", 
            rho=as.double(rho),
            M=as.integer(M),
            G = as.integer(G),
            P = as.integer(P),
            g = as.integer(g),
            p = as.integer(p),
            theta = as.double(theta),
            mu = as.double(mu),
            sumtheta = as.double(0)
            )
  return(out)
}