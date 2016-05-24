######  R package

### from data compute the MLEs in the exponential case
.MLE_exp = function(ydat,xdat) {
              n = length(ydat)
              m = length(xdat)
	        my = mean(xdat)
	        mx = mean(ydat)
              l_hat = (1/my) + (1/mx)
              psi_hat = my/(mx+my) 
              return(cbind(l_hat,psi_hat))
              }


### from data compute the MLEs in the nornal case with equal variances
.MLE_normEV = function(ydat,xdat) {
              n = length(ydat)
              m = length(xdat)
	        my = mean(xdat)
	        mx = mean(ydat)
              sigma1_hat = (1/n)*(sum((ydat-mx)^2))
              sigma2_hat = (1/m)*(sum((xdat-my)^2))
              sigma_hat = (n*sigma1_hat + m*sigma2_hat)/(n+m)
              lambda1_hat = mx/sqrt(2*sigma_hat)
              lambda2_hat = sqrt(2*sigma_hat)
              psi_hat = pnorm((my-mx)/sqrt(2*sigma_hat))
              Di_hat = lambda2_hat * (qnorm(psi_hat) + lambda1_hat)
              return(cbind(lambda1_hat,lambda2_hat,psi_hat,Di_hat))
           }



### from data compute the MLEs in the nornal case with different variances
.MLE_normDV = function(ydat,xdat) {
              n = length(ydat)
              m = length(xdat)
	        my = mean(xdat)
	        mx = mean(ydat)
              sigma1_hat = (1/n)*(sum((ydat-mx)^2))
              sigma2_hat = (1/m)*(sum((xdat-my)^2))
              lambda1_hat = mx
              lambda2_hat = sigma1_hat
              lambda3_hat = sigma2_hat
              psi_hat = pnorm((my-mx)/sqrt(sigma1_hat+sigma2_hat))
              b_hat = sqrt(lambda2_hat + lambda3_hat)
              Di_hat = b_hat*qnorm(psi_hat) + lambda1_hat
              return(cbind(lambda1_hat,lambda2_hat,lambda3_hat,psi_hat,Di_hat,b_hat))
           }




MLEs = function(ydat,xdat,distr) {
       switch(distr, 
             exp = .MLE_exp(ydat,xdat),
             norm_EV = .MLE_normEV(ydat,xdat),
             norm_DV = .MLE_normDV(ydat,xdat))
       }




loglik =  function(ydat,xdat,lambda,psi,distr="exp"){
          n = length(ydat)
          m = length(xdat)
          if (distr!="exp") {
                if (distr=="norm_EV") {
                      MLE=.MLE_normEV(ydat,xdat)
                      l1_hat <- MLE[1]
                      l2_hat <- MLE[2]
                      psi_hat <- MLE[3]
                      Di_hat <- MLE[4]
                      Di = lambda[2] * (qnorm(psi) + lambda[1])
                      loglikelihood=-(n+m)*log(lambda[2]) - (1/lambda[2]^2)*( (n+m)*l2_hat^2/2 +
                             n*(lambda[1]*lambda[2] - l1_hat*l2_hat)^2 + m*(Di-Di_hat)^2 )                             
                } else { 
                      if (distr=="norm_DV") {
                             MLE=.MLE_normDV(ydat,xdat)
                             l1_hat <- MLE[1]
                             l2_hat <- MLE[2]
                             l3_hat <- MLE[3]
                             psi_hat <- MLE[4]
                             Di_hat <- MLE[5]
                             b_hat <- MLE[6]
                             Di = (sqrt(lambda[2] + lambda[3]))*qnorm(psi) + lambda[1] 
                             qu=(lambda[1]-l1_hat)^2
                             loglikelihood= -(1/2)*(n*log(lambda[2])+m*log(lambda[3])) - 
                                            (n/(2*lambda[2]))*(l2_hat + qu)-(m/(2*lambda[3]))*(l3_hat +(Di_hat-Di)^2)                                     
                       } else {
                             stop("Error: distribution not valid, choose between 'exp', 'norm_ev' and 'norm_DV'")
                              }
                       }
             } else {
                  my = mean(xdat)
	            mx = mean(ydat)
                  psi_hat=my/(mx+my)
                  l_hat= 1/mx + 1/my 
                  MLE= c(l_hat,psi_hat)
                  loglikelihood= n*log(psi)+m*log(1-psi)+(n+m)*log(lambda)-
                                 n*psi*lambda*mx - m*(1-psi)*lambda*my
             }
             return(list(psi_hat=psi_hat,log_likelihood =loglikelihood))
             }





.loglik2 =  function(N,MLE,lambda,psi,distr="exp"){
           n=N[1]
           m=N[2]
           if (distr!="exp") {
                if (distr=="norm_EV") {
                      l1_hat <- MLE[1]
                      l2_hat <- MLE[2]
                      psi_hat <- MLE[3]
                      Di_hat <- MLE[4]
                      Di = lambda[2] * (qnorm(psi) + lambda[1])
                      loglikelihood=-(n+m)*log(lambda[2]) - (1/lambda[2]^2)*( (n+m)*l2_hat^2/2 +
                                    n*(lambda[1]*lambda[2] - l1_hat*l2_hat)^2 + m*(Di-Di_hat)^2 )                             
                } else { 
                      if (distr=="norm_DV") {
                             l1_hat <- MLE[1]
                             l2_hat <- MLE[2]
                             l3_hat <- MLE[3]
                             psi_hat <- MLE[4]
                             Di_hat <- MLE[5]
                             b_hat <- MLE[6]
                             Di = (sqrt(lambda[2] + lambda[3]))*qnorm(psi) + lambda[1] 
                             qu=(lambda[1]-l1_hat)^2
                             loglikelihood= -(1/2)*(n*log(lambda[2])+m*log(lambda[3])) -
                                    (n/(2*lambda[2]))*(l2_hat + qu)-(m/(2*lambda[3]))*(l3_hat +(Di_hat-Di)^2) 
                      } else {
                            stop("Error: distribution not valid, choose between 'exp', 'norm_ev' and 'norm_DV'")
                              }
                       }
             } else {
                  l_hat= MLE[1] 
                  psi_hat=MLE[2]
                  loglikelihood= n*log(psi)+m*log(1-psi)+(n+m)*log(lambda)-
                         n*psi*lambda/(l_hat*psi_hat) - m*(1-psi)*lambda/(l_hat*(1-psi_hat))
             }
             return(list(psi_hat=psi_hat,log_likelihood =loglikelihood))
            }





rp =function(ydat,xdat,psi,distr="exp") {
                results = .nuisance(ydat,xdat,psi,distr)
                N = results[[1]]
                mle = results[[2]]
                lambda_psi= results[[3]]
                lambda_hat= mle[1:(length(mle)/2)]
                psi_hat = mle[(length(mle)/2)+1]
                p_loglik= .loglik2(N,mle,lambda_psi,psi,distr)[[2]]    # profile log-likelihood for given psi
                p_loglik_psihat = .loglik2(N,mle,lambda_hat,psi_hat,distr)[[2]]   # log-likelihood in the MLE                
                return(list( 
                   signed_loglikelihood_ratio = sign(psi_hat-psi)*sqrt(2*(p_loglik_psihat - p_loglik))
                    )) 
     }




.rp2 =function(ydat,xdat,psi,distr="exp") {
                results = .nuisance(ydat,xdat,psi,distr)
                N = results[[1]]
                mle = results[[2]]
                lambda_psi= results[[3]]
                lambda_hat= mle[1:(length(mle)/2)]
                psi_hat = mle[(length(mle)/2)+1]
                p_loglik= .loglik2(N,mle,lambda_psi,psi,distr)[[2]]    # profile log-likelihood for given psi
                p_loglik_psihat = .loglik2(N,mle,lambda_hat,psi_hat,distr)[[2]]   # log-likelihood in the MLE                
                return(list( 
                   N=N, MLE=mle, Lambda_psi=lambda_psi,
                   signed_loglikelihood_ratio = sign(psi_hat-psi)*sqrt(2*(p_loglik_psihat - p_loglik))
                    )) 
     }




.nuisance = function(ydat,xdat,psi,distr="exp") {
           n=length(ydat)
           m=length(xdat)
           if (distr!="exp") {
                if (distr=="norm_EV") {
                       mle=.MLE_normEV(ydat,xdat)
                       # c("l1_hat","l2_hat","psi_hat","Di_hat")
                       log_lik =function(lambda){
                                 Di = lambda[2] * (qnorm(psi) + lambda[1])
                                 -( -(n+m)*log(lambda[2]) - (1/lambda[2]^2)*( (n+m)*mle[2]^2/2 +
                                 n*(lambda[1]*lambda[2] - mle[1]*mle[2])^2 + m*(Di-mle[4])^2 ) )
                                }
                       O<-optim( c(3,2), log_lik, NULL, method="L-BFGS-B", lower=c(-Inf,0.0001), hessian = TRUE)    
                       lam_psi=O$par
                } else { 
                      if (distr=="norm_DV") {
                             mle=.MLE_normDV(ydat,xdat)
                             # c("l1_hat","l2_hat","l3_hat","psi_hat","Di_hat","b_hat")
                             log_lik = function(lambda){
                                       Du = (sqrt(lambda[2] + lambda[3]))*qnorm(psi) + lambda[1] 
                                       qu=(lambda[1]-mle[1])^2
                                       -(1/2)*(n*log(lambda[2])+m*log(lambda[3])) -(n/(2*lambda[2]))*(mle[2] + qu)-
                                       (m/(2*lambda[3]))*(mle[3] +(mle[5]-Du)^2) 
                                      }
                             log_lik_g = function(lambda){
                                          if ((lambda[2]<=0.001)|(lambda[3]<=0.001)) { -1000000000000
                                          } else {log_lik(lambda)}  }
                             O<-optim( c(5,1,1.5), log_lik_g, NULL, method = "BFGS", control = list(fnscale=-1),hessian = TRUE)
                             lam_psi=O$par
                       } else {
                            stop("Error: distribution not valid, choose between 'exp', 'norm_ev' and 'norm_DV'")
                              }
                 }
           } else {
           mle =.MLE_exp(ydat,xdat)
           lam_psi = ((n+m)*mle[1])/((n*psi/mle[2])+(m*(1-psi)/(1-mle[2]))) 
                   }
           return(list(N=c(n,m), MLE=mle, lambda_psi=lam_psi))
          }





rpstar = function(ydat,xdat,psi,distr="exp") {
                results =.rp2(ydat,xdat,psi,distr) 
                n = results[[1]][1]
                m = results[[1]][2]
                mle = results[[2]]
                lambda_psi= results[[3]]
                r_p = results[[4]]
                if (distr!="exp") {
                     if (distr=="norm_EV") {
                           l1_hat <- mle[1]
                           l2_hat <- mle[2]
                           psi_hat <- mle[3]
                           Di_hat <- mle[4]
                           # l_tilde /(J_hat *J_lambda_psi)^1/2
                           b = lambda_psi[2]*(qnorm(psi) + lambda_psi[1])
                           d_psihat = 1/dnorm(qnorm(psi_hat))
                           num =(2*l2_hat^2)*(2*n*m*(l1_hat*l2_hat-Di_hat)^2 + ((n+m)*l2_hat)^2 -
                                 n*m*(l1_hat*l2_hat-Di_hat)*(lambda_psi[1]*lambda_psi[2]-b))
                           den = lambda_psi[2]^4*(12*n*m*(l1_hat*l2_hat-Di_hat)^2 + 
                                 8*(n*l1_hat*l2_hat + m*Di_hat)^2 - 2*(n+m)^2*(lambda_psi[2]^2-3*l2_hat^2) -
                                 8*(n+m)*(n*lambda_psi[1]*lambda_psi[2]*l1_hat*l2_hat +m*Di_hat*b))
                           A= 2*n*m*(n+m)*l2_hat^2 * d_psihat^2
                           q2 = num/sqrt(den*A)
                           # |num| in q(psi)
                           q1 = .lderiv_nEV(n,m,l1_hat,l2_hat,Di_hat,d_psihat,lambda_psi,b)  # matrices of drivatives of log-lik
                           qu= q1 * q2                                                                             
                      } else {
                           if (distr=="norm_DV") {
                                 l1_hat <- mle[1]
                                 l2_hat <- mle[2]
                                 l3_hat <- mle[3]
                                 psi_hat <- mle[4]
                                 Di_hat <- mle[5]
                                 b_hat <- mle[6]                
                                 qu = .lderiv_nDV(n,m,psi,l1_hat,l2_hat,l3_hat,psi_hat,Di_hat,b_hat,lambda_psi)                                   
                           } else {
                              stop("Error: distribution not valid, choose between 'exp', 'norm_ev' and 'norm_DV'")
                                  }
                             }
                  } else {
                      qu = .lderiv_exp(n,m,psi,mle)   
                         }
                  if (is.na(qu)) {warning("Warning: could not compute statistics because of imprecise estimates. Increase sample sizes slightly.")}
                  rp_s =r_p+(1/r_p)*log(abs(qu/r_p))     
                  return(list(rp= r_p, rp_star= rp_s))
          }





.lderiv_nEV = function(n,m,l1_hat,l2_hat,Di_hat,d_psihat,lambda_psi,b) {
                     # ltilde_lambda;lambdahat
                     lt_11 = 2*(n+m)* l2_hat /lambda_psi[2]
                     lt_12 = (2/lambda_psi[2])* (n*l1_hat + m*Di_hat/l2_hat)
                     lt_21 = ((2*l2_hat/(lambda_psi[2]^3))* (2*n*l1_hat*l2_hat +
                             2*m*Di_hat - m*b - n*lambda_psi[1]*lambda_psi[2] ) )
                     lt_22 = ( (2/(lambda_psi[2]^3))* ((n+m)*l2_hat + 2*n*(l1_hat^2)*l2_hat +
                             2*m*Di_hat^2/l2_hat - n*lambda_psi[1]*lambda_psi[2]*l1_hat - m*b*Di_hat/l2_hat) )
                     ltilde = matrix(c(lt_11,lt_21,lt_12,lt_22),nrow = 2, ncol = 2)
                     # lhat_ ;lambdahat  -  ltilde_ ;lambdahat
                     v3 = matrix(c(0, -(n+m)/l2_hat),nrow = 2, ncol = 1)
                     v21 = -2*(l2_hat/lambda_psi[2]^2)*(n*(l1_hat*l2_hat-lambda_psi[1]*lambda_psi[2])+m*(Di_hat-b))
	               v22 =( -(1/lambda_psi[2]^2)*((n+m)*l2_hat+2*n*(l1_hat^2)*l2_hat +
                          2*m*(Di_hat-b)*Di_hat/l2_hat -2*n*l1_hat*lambda_psi[1]*lambda_psi[2])  )
	               v2 = matrix(c(v21,v22), nrow = 2, ncol = 1)
                     # ltilde_ lambda;psihat
                     v11 = 2*m*l2_hat*d_psihat/lambda_psi[2]
	               v12 = (2*m*l2_hat*d_psihat/lambda_psi[2]^3)*(2*Di_hat - b) 
	               v1 = matrix(c(v11,v12), nrow = 1, ncol = 2)
                     # ltilde_ lambda;psihat
                     n2 = -2*m*(Di_hat-b)*l2_hat*d_psihat/lambda_psi[2]^2
                     # numer of the q term
                     num = -n2-v1%*%solve(ltilde)%*%(v3-v2)

                     return( num )
             }

 




.lderiv_nDV = function(n,m,psi,l1_hat,l2_hat,l3_hat,psi_hat,Di_hat,b_hat,lambda_psi) {
                     b_til = sqrt(lambda_psi[2] + lambda_psi[3])
                     Di_til = b_til*qnorm(psi) + lambda_psi[1]
                     d_psihat = 1/dnorm(qnorm(psi_hat))  # first derivative of qnorm(psi_hat)
                     # ltilde_lambda;lambdahat
                     lt_11 = n/lambda_psi[2] + m/lambda_psi[3]
                     lt_12 = qnorm(psi_hat)*m/(2*lambda_psi[3]*b_hat)
                     lt_13 = qnorm(psi_hat)*m/(2*lambda_psi[3]*b_hat)
                     lt_21 = (l1_hat-lambda_psi[1])*n/lambda_psi[2]^2 + m*qnorm(psi)/(2*lambda_psi[3]*b_til)
                     lt_22 = (1/2)*(n/lambda_psi[2]^2 + m*qnorm(psi)*qnorm(psi_hat)/(2*lambda_psi[3]*b_hat*b_til))
                     lt_23 = m*qnorm(psi)*qnorm(psi_hat)/(4*lambda_psi[3]*b_hat*b_til)
                     lt_31 = m*(Di_hat-Di_til)/lambda_psi[3]^2 + m*qnorm(psi)/(2*lambda_psi[3]*b_til)
                     lt_32 = m*(Di_hat-Di_til)*qnorm(psi_hat)/(2*lambda_psi[3]^2*b_hat) + m*qnorm(psi)*qnorm(psi_hat)/(4*lambda_psi[3]*b_hat*b_til)
                     lt_33 = m/(2*lambda_psi[3]^2)*(1 + qnorm(psi_hat)*(Di_hat-Di_til)/b_hat) + m*qnorm(psi)*qnorm(psi_hat)/(4*lambda_psi[3]*b_hat*b_til)
                     l_tilde = matrix(c(lt_11,lt_12,lt_13,lt_21,lt_22,lt_23,lt_31,lt_32,lt_33),nrow=3, ncol=3,byrow=TRUE)
                     # J_lambda,lambda
                     # p=vector of parameters  p[1]= lambda1, p[2]= lambda2, p[3]= lambda3, p[4] = psi
                     J_compute = function(p){
                               bb = sqrt(p[2] + p[3])
                               DD = bb*qnorm(p[4]) + p[1]
                               l_11 = -n/p[2] - m/p[3]
                               l_12 = -(l1_hat-p[1])*n/(p[2])^2 - m*qnorm(p[4])/(2*p[3]*bb)
                               l_13 = -m*qnorm(p[4])/(2*p[3]*bb) - m*(Di_hat-DD)/(p[3])^2
                               l_22 = n/(2*(p[2])^2) - n*(l2_hat + (p[1]-l1_hat)^2)/(p[2])^3 + m*qnorm(p[4])*(DD-Di_hat)/(4*p[3]*bb^3) - m*(qnorm(p[4]))^2/(4*p[3]*bb^2)
                               l_23 = - m*(qnorm(p[4]))^2/(4*p[3]*bb^2) + m*qnorm(p[4])*(2*p[2]+3*p[3])*(DD-Di_hat)/(4*(p[3])^2*bb^3)
                               l_33 = m*(1/2 - l3_hat/p[3] - (Di_hat-DD)^2/p[3])/(p[3])^2 - m*(qnorm(p[4]))^2/(4*p[3]*bb^2) - 
                               m*qnorm(p[4])*(4*p[2]+5*p[3])*(Di_hat-DD)/(4*(p[3])^2*bb^3)
                               matrix(c(-l_11,-l_12,-l_13,-l_12,-l_22,-l_23,-l_13,-l_23,-l_33),nrow=3,ncol=3,byrow=TRUE)
                              }
                     # Jhat_lambda,lambda
                     Jll_hat =J_compute(c(l1_hat,l2_hat,l3_hat,psi_hat))  
                     # Jtilde_lambda,lambda
                     Jll_tilde =J_compute(c(lambda_psi,psi))  
                     # l_tilde ;psi_hat
                     n2 = m*b_hat*d_psihat*(Di_til-Di_hat)/lambda_psi[3]
                     # l_tilde lambda;psi_hat
                     v11 = m*b_hat*d_psihat/lambda_psi[3]
                     v12 = m*qnorm(psi)*b_hat*d_psihat/(2*lambda_psi[3]*b_til)
                     v13 = m*b_hat*d_psihat*(Di_hat-Di_til)/(lambda_psi[3]^2) + m*qnorm(psi)*b_hat*d_psihat/(2*lambda_psi[3]*b_til)
                     v1 = matrix(c(v11,v12,v13), nrow = 1, ncol = 3)
                     # l_hat ;lambda_hat
                     l_hat1 = 0
                     l_hat2 = -n/(2*l2_hat) 
                     l_hat3 = -m/(2*l3_hat)  
                     v3= matrix(c(l_hat1,l_hat2,l_hat3),nrow = 3, ncol = 1)
                     # l_tilde ;lambda_hat
                     v21 = n*(lambda_psi[1]-l1_hat)/lambda_psi[2] - m*(Di_hat-Di_til)/lambda_psi[3]
                     v22 = -n/(2*lambda_psi[2]) - m*qnorm(psi_hat)*(Di_hat-Di_til)/(2*lambda_psi[3]*b_hat)
                     v23 = -m*qnorm(psi_hat)*(Di_hat-Di_til)/(2*lambda_psi[3]*b_hat) - m/(2*lambda_psi[3])
                     v2 =  matrix(c(v21,v22,v23), nrow = 3, ncol = 1)
                     # J_hat_psi,psi
                     jp <- m*(b_hat^2)*(d_psihat^2)/l3_hat
                     # J_hat_psi,lambda
                     jpl1 = m*b_hat*d_psihat/l3_hat
                     jpl2 = m*qnorm(psi_hat)*d_psihat/(2*l3_hat)
                     jpl3 = m*qnorm(psi_hat)*d_psihat/(2*l3_hat)
                     jpl = matrix(c(jpl1,jpl2,jpl3),nrow=1,ncol=3)
                     q1 = -n2-v1%*%solve(l_tilde)%*%(v3-v2)
	               den = jp-jpl%*%solve(Jll_hat)%*%t(jpl)   
                    # if ((det(Jll_tilde)*det(Jll_hat)*den)<0) { return(NA)}
                     den_q2 <- det(Jll_tilde)*det(Jll_hat)*den
                   # if ((det(Jll_tilde)*det(Jll_hat)*den)<0) { return(NA)}
                    if (den_q2 <0 ) den_q2 <- abs( det(Jll_tilde)*det(Jll_hat)*den )
                     q2 = det(l_tilde)/sqrt(den_q2)
                     return(q1*q2)
            }
       



.lderiv_exp = function(n,m,psi,mle) {     
                      l_hat= mle[1] 
                      psi_hat=mle[2]
                      l_psi =((n+m)*l_hat)/((n*psi/psi_hat)+(m*(1-psi)/(1-psi_hat)))   ## lambda hat given psi
                      #l;psihat ,  in p=vettore dei parametri p[1]= psi, p[2]= lambda
                      sample_der1 = function(p){
                                  ((n*p[1]*p[2])/(psi_hat^2*l_hat)) - ((m*p[2]*(1-p[1]))/(l_hat*(1-psi_hat)^2)) }
                      inp=list(c(psi_hat,l_hat),c(psi,l_psi))
                      l_mle =unlist(lapply(inp,sample_der1))  #l;psihat on psi_hat and on l_psi
                      #l_lambda;mle
                      l_21 =(1/l_hat)*((n*psi/psi_hat^2)-(m*(1-psi)/((1-psi_hat)^2)))
                      #l_lambda;lambdahat
                      l_22 =(1/l_hat^2)*((n*psi/psi_hat)+(m*(1-psi)/(1-psi_hat)))
                      #l;lambdahat
                      sample_der2 = function(p){
                                   ((n*p[1]*p[2])/(psi_hat*l_hat^2)) + ((m*p[2]*(1-p[1]))/(l_hat^2*(1-psi_hat))) }
                      l_lamhat =unlist(lapply(inp,sample_der2))  #l;lambdahat on psihat and on l_psi
                      # det(J(thetahat))
                      J= n*m/(l_hat^2*psi_hat^2*((1-psi_hat)^2))
                      #j_lambda_lambda
                      j22 =(n+m)/l_psi^2
                      ## q
                      qu = (l_mle[1]-l_mle[2]-((l_21/l_22)*(l_lamhat[1]-l_lamhat[2])))*l_22/(sqrt(J*j22))
                      return(qu)
              }




wald = function(ydat,xdat,psi,distr="exp") { 
                results = .nuisance(ydat,xdat,psi,distr)
                n = results[[1]][1]
                m = results[[1]][2]
                mle = results[[2]]
                lambda_psi= results[[3]]
                if (distr!="exp") {
                     if (distr=="norm_EV") {
                           l1_hat <- mle[1]
                           l2_hat <- mle[2]
                           psi_hat <- mle[3]
                           Di_hat <- mle[4]
                           # first and second derivatives of qnorm(psi_hat)
                           d_psihat = 1/dnorm(qnorm(psi_hat))
                           d2_psihat = (qnorm(psi_hat)) * (d_psihat^2)
                           # jhat_lambda
                           C = ((n+m)/2)* (l2_hat)^2 + n* (l1_hat)^2 * (l2_hat)^2 + m*Di_hat^2
                           jh_11 = 2*(n+m)
                           jh_22 = - (n+m+4*n*(l1_hat^2))/(l2_hat^2) - (4*m*Di_hat^2)/(l2_hat^4) + 6*C/(l2_hat^4)
                           jh_12 =  2*n*l1_hat/l2_hat + 2*m*Di_hat/(l2_hat^2)
                           j_hat = matrix(c(jh_11,jh_12,jh_12,jh_22),nrow = 2, ncol = 2)
                           # J_hat_psi,psi
                           jpp = 2*m*(d_psihat)^2
                           # J_hat_psi,lambda                        
	                     jpl1 = 2*m*d_psihat
	                     jpl2 = (2*m/l2_hat)*Di_hat*d_psihat/l2_hat
                           jpl = matrix(c(jpl1,jpl2), nrow = 1, ncol = 2)
                           # Wald statistic
                           jp_hat = abs(jpp-jpl%*%solve(j_hat)%*%t(jpl))
                           w = (psi_hat-psi)*sqrt(jp_hat)                                                                                              
                      } else {
                           if (distr=="norm_DV") {
                                 l1_hat <- mle[1]
                                 l2_hat <- mle[2]
                                 l3_hat <- mle[3]
                                 psi_hat <- mle[4]
                                 Di_hat <- mle[5]
                                 b_hat <- mle[6]      
                                 # first  derivative of qnorm(psi_hat)   
                                 d_psihat = 1/dnorm(qnorm(psi_hat))                      
                                 # J_hat_psi,psi
                                 jp <- m*(b_hat^2)*(d_psihat^2)/l3_hat
                                 # J_hat_psi,lambda
                                 jpl1 = m*b_hat*d_psihat/l3_hat
                                 jpl2 = m*qnorm(psi_hat)*d_psihat/(2*l3_hat)
                                 jpl3 = m*qnorm(psi_hat)*d_psihat/(2*l3_hat)
                                 jpl = matrix(c(jpl1,jpl2,jpl3),nrow=1,ncol=3)
                                 # p=vector of parameters  p[1]= lambda1, p[2]= lambda2, p[3]= lambda3, p[4] = psi
                                 J_compute1 = function(p){
                                             bb = sqrt(p[2] + p[3])
                                             DD = bb*qnorm(p[4]) + p[1]
                                             l_11 = -n/p[2] - m/p[3]
                                             l_12 = -(l1_hat-p[1])*n/(p[2])^2 - m*qnorm(p[4])/(2*p[3]*bb)
                                             l_13 = -m*qnorm(p[4])/(2*p[3]*bb) - m*(Di_hat-DD)/(p[3])^2
                                             l_22 = n/(2*(p[2])^2) - n*(l2_hat + (p[1]-l1_hat)^2)/(p[2])^3 + m*qnorm(p[4])*(DD-Di_hat)/(4*p[3]*bb^3) - m*(qnorm(p[4]))^2/(4*p[3]*bb^2)
                                             l_23 = - m*(qnorm(p[4]))^2/(4*p[3]*bb^2) + m*qnorm(p[4])*(2*p[2]+3*p[3])*(DD-Di_hat)/(4*(p[3])^2*bb^3)
                                             l_33 = m*(1/2 - l3_hat/p[3] - (Di_hat-DD)^2/p[3])/(p[3])^2 - m*(qnorm(p[4]))^2/(4*p[3]*bb^2) - 
                                                    m*qnorm(p[4])*(4*p[2]+5*p[3])*(Di_hat-DD)/(4*(p[3])^2*bb^3)
                                             matrix(c(-l_11,-l_12,-l_13,-l_12,-l_22,-l_23,-l_13,-l_23,-l_33),nrow=3,ncol=3,byrow=TRUE)
                                            }
                                 # Jhat_lambda,lambda
                                 Jll_hat =J_compute1(c(l1_hat,l2_hat,l3_hat,psi_hat))  
                                 # Wald statistic
                                 jp_hat = jp-jpl%*%solve(Jll_hat)%*%t(jpl)
                                 w = (psi_hat-psi)*sqrt(jp_hat)
                           } else {
                              stop("Error: distribution not valid, choose between 'exp', 'norm_EV' and 'norm_DV'")
                                  }
                             }
                  } else {
                      l_hat= mle[1] 
                      psi_hat=mle[2]
                      # Wald statistic
                      jp_hat = (n*m)/((n+m)*(psi_hat^2*(1-psi_hat)^2))
                      w = (psi_hat-psi)*sqrt(jp_hat)
                         }
                  return(list(Wald=w, Jp_hat=jp_hat))
          }





Prob = function(ydat,xdat,distr="exp", method="RPstar",level=0.05){
                z = qnorm(level/2)
                if (method !="Wald") { 
                      pp =seq(0.001,0.999,by=0.001)  ## values for psi
                      myfun=rep(0,length(pp))
                      if (method=="RP") {
                            for (i in 1:length(pp)) {
                                  myfun[i] = rp(ydat,xdat,pp[i],distr)[[1]]  }                          
                      } else { if (method=="RPstar") {
                            for (i in 1:length(pp)) {
                                  myfun[i]= rpstar(ydat,xdat,pp[i],distr)[[2]]   }
                                     # myfun = sapply(pp,rpstar,ydat,xdat,distr)     # maybe, in rpstar pp should be the first element passed
                                } else { stop("Error: method is not valid, choose 'Wald', 'RP' or 'RPstar'")  }
                      } 
                      rsmooth = smooth.spline(myfun,pp)
                      # point estimate norm_EV"
                      rprev = predict(rsmooth,x=0)
                      psi_hat = rprev$y
                      # Confidence interval
                      UB = predict(rsmooth,x=z)$y 
                      LB = predict(rsmooth,x=-z)$y
                } else { 
                      jp_hat = wald(ydat,xdat,0.5,distr)[[2]]  
                       if (distr=="exp"){                      
                             psi_hat = .MLE_exp(ydat,xdat)[2]
                       } else { if (distr=="norm_EV") {
                                     psi_hat = .MLE_normEV(ydat,xdat)[3]
                                } else {
                                     psi_hat = .MLE_normDV(ydat,xdat)[4]
                                }
                       }
                      UB = psi_hat-(jp_hat^(-0.5))*z
                      LB = psi_hat+(jp_hat^(-0.5))*z
                }
                return(list(PROB=psi_hat, C.Interval=c(LB,UB) ))
        }
              



ROC.plot = function(ydat,xdat,distr="exp", method="RPstar",mc=1){
                    if ((method == "RPstar")|(mc!=1)){
                         a = Prob(ydat,xdat,distr,"RPstar")
                         res = .nuisance(ydat,xdat,a[[1]],distr)    ## with psi=psi_hat = a[[1]],  if method=="RPstar" then a[[1]]=psi_hat^* 
                         lamhat_star =res[[3]]        
                         if (distr=="exp"){
                             alpha = a[[1]]*lamhat_star         ## psi_hat^* * lambda_hat^*
                             beta = lamhat_star*(1-a[[1]])
                             t1=seq(0,50,by=0.01)  ## values for the cut-off 
                             sens_1=pexp(t1,rate=beta,lower.tail=FALSE) 
                             spec_1 =pexp(t1,rate=alpha,lower.tail=FALSE)
                         } else {
                             if (distr=="norm_EV"){
                                 mux = lamhat_star[1] * lamhat_star[2]    ## mux= lambda1_hat^* * lambda2_hat^* 
                                 muy = lamhat_star[2]*(qnorm(a[[1]]) + lamhat_star[1]) 
                                 sigma2x = sigma2y = (lamhat_star[2]^2)/2
                             } else {
                                 mux = lamhat_star[1]
                                 muy = (sqrt(lamhat_star[2]+lamhat_star[3]))*qnorm(a[[1]]) + lamhat_star[1]
                                 sigma2x = lamhat_star[2]
                                 sigma2y = lamhat_star[3]
                             }
                             range=c(min(mux,muy)-5*max(sqrt(sigma2x),sqrt(sigma2y)), max(mux,muy)+5*max(sqrt(sigma2x),sqrt(sigma2y)))
                             t1=seq(range[1],range[2],by=0.01)  ## values for the cut-off 
                             sens_1=pnorm(t1,muy,sqrt(sigma2y),lower.tail=FALSE)
                             spec_1 =pnorm(t1,mux,sqrt(sigma2x),lower.tail=FALSE)
                         }
                      }
                      if ((method!="RPstar")|(mc!=1)) {
                         mle = MLEs(ydat,xdat,distr)
                         if (distr=="exp"){
                             alpha = mle[2]*mle[1]           ## psi_hat * lambda_hat
                             beta = mle[1]*(1-mle[2])
                             t2=seq(0,50,by=0.01)  ## values for the cut-off 
                             sens_2=pexp(t2,rate=beta,lower.tail=FALSE)
                             spec_2 =pexp(t2,rate=alpha,lower.tail=FALSE)
                         } else {
                             if (distr=="norm_EV"){
                                 mux = mle[1]* mle[2]    ## mux= lambda1_hat * lambda2_hat 
                                 muy = mle[4] 
                                 sigma2x = sigma2y = (mle[2]^2)/2
                             } else {
                                 mux = mle[1]
                                 muy = mle[5]
                                 sigma2x = mle[2]
                                 sigma2y = mle[3]
                             }
                             range=c(min(mux,muy)-5*max(sqrt(sigma2x),sqrt(sigma2y)), max(mux,muy)+5*max(sqrt(sigma2x),sqrt(sigma2y)))
                             t2=seq(range[1],range[2],by=0.01)  ## values for the cut-off 
                             sens_2=pnorm(t2,muy,sqrt(sigma2y),lower.tail=FALSE)
                             spec_2 =pnorm(t2,mux,sqrt(sigma2x),lower.tail=FALSE)
                         }
                    }
                    plot(-1,-1,type="l",lwd=2,xlim=c(0,1), ylim=c(0,1),xlab="1- specificity",ylab="sensitivity")
                    title("ROC curve")
                    lines(c(0,1),c(0,1))
                    if (mc==1){
                        if (method=="RPstar"){
                             lines(spec_1,sens_1,type="l",lwd=2)
                        } else { 
                             lines(spec_2,sens_2,type="l",lwd=2)
                        }
                    } else {
                        lines(spec_1,sens_1,type="l",lwd=2,col=1)
                        lines(spec_2,sens_2,type="l",lwd=2,col=2)
                        legend(0.7, 0.2, c("r*","MLE"), col = c(1,2),lwd=c(2,2))
                    }
             }

                          
                             
                    
                        
                         
                     
