# function to estimate jointly the attachment function and node fitness  2015-3-11 Thong Pham
PAFit <- function(data, only_PA = FALSE, only_f = FALSE, mode_f = c("Linear_PA", "Constant_PA"),
                  true_A = NULL, true_f = NULL,
                  shape = 0.1 , rate = 0.1, 
                  auto_lambda = TRUE, ratio = 0, 
                  lambda = 1, weight_PA_mode = c(0,1), auto_stop = TRUE,stop_cond = 10^-4, iteration = 20, 
                  max_iter = 1000,
                  debug = FALSE,
                  step_size = 1, q = 1, 
                  normalized_f = FALSE, interpolate = FALSE,...) {
    if ((data$only_PA == TRUE) & (only_PA == FALSE)) {
        stop("Error: the data do not support estimation of node fitness. Please re-run GetStatistics again with the option 'only_PA = FALSE', or run PAFit with option 'only_PA = TRUE'")
          
    }
    if ((shape <= 0) || (rate <= 0))
        stop("Error: shape and rate should be positive number")
  
    non_zero_theta     <- which(data$Sum_m_k > 0)
    num_nonzero        <- length(non_zero_theta)
    theta              <- rep(1,length(data$Sum_m_k))
    
    
    
     for (ii in 1:length(theta))
        if ((mode_f[1] == "Constant_PA") & (only_f == TRUE)) 
        {
            theta[ii] <- 1;  
        } else if ((only_f == TRUE) & (!is.null(true_A))) {
            theta[ii] <- true_A[ii];    
        }
        else {    
            if (ii <= data$start_deg & ii %in% non_zero_theta)
               theta[ii] <- ii - 1
            else 
            if (ii > data$start_deg & ii %in% non_zero_theta)
               theta[ii] <- data$begin_deg[ii - data$start_deg] - 1
        }
          
    theta[which(theta <= 0)] <- 1
    theta[-non_zero_theta] <- 0
    theta         <- theta/sum(theta)
    non_zero_f    <- which(data$z_j != 0)
    f             <- rep(2,length(data$f_position))
    f             <- length(f) * f/sum(f)
    
    if (TRUE == only_PA) {
        if (is.null(true_f))
            f[] <- 1
        else
            f[] <- true_f  
    }
    log_likelihood   <- vector()
    
    update_theta   <- non_zero_theta[-c(1,2,num_nonzero,num_nonzero - 1)]
    noupdate_theta <- non_zero_theta[c(1,2,num_nonzero,num_nonzero - 1)]
    minus_1 <- non_zero_theta[-c(1,num_nonzero,num_nonzero - 1,num_nonzero-2)]
    minus_2 <- non_zero_theta[-c(num_nonzero,num_nonzero - 1,num_nonzero-2,num_nonzero-3)] 
    plus_1  <- non_zero_theta[-c(1,2,3,num_nonzero)]
    plus_2  <- non_zero_theta[-c(1,2,3,4)] 

    
    #weights of the regularization term
    if (weight_PA_mode[1] == 0)
        w_k <- data$Sum_m_k/sum(data$Sum_m_k)
    else
        w_k <- rep(1/length(theta),length(theta))
    
    if (TRUE == auto_lambda) {
      lambda <- ratio * sum(data$Sum_m_k)
    }
    if (TRUE == auto_stop)
        iteration <- max_iter
    if (q > 1) {
       parameter_save <- matrix(0,nrow = num_nonzero + length(non_zero_f), ncol = q)
       U              <- matrix(0,nrow = num_nonzero + length(non_zero_f), ncol = q - 1)
       V              <- matrix(0,nrow = num_nonzero + length(non_zero_f), ncol = q - 1)
       candidate      <- rep(0,num_nonzero + length(non_zero_f))
    }
    
    candidate_ok     <- 0
    candidate_accept <- 0
    normalized_const <- rep(0, dim(data$n_tk)[1])
    for (i in 1:iteration) {
      
        if ((!is.null(true_f)) || (FALSE == only_PA)){
            .normalized_constant(normalized_const,data$node_degree,theta,f,data$offset_tk) 
        } else if (is.null(true_f)) {
            normalized_const <- as.vector(data$n_tk[,non_zero_theta]%*%theta[non_zero_theta])
        } else stop("Not yet implemented functional")
        
        time_non_zero     <- which(normalized_const != 0)
        non_zero_f_temp   <- which(f >= 0)
        if ((FALSE == only_PA) && (FALSE == only_f)) {
            log_likelihood    <- c(log_likelihood, sum(data$z_j[non_zero_f_temp] * log(f[non_zero_f_temp])) + 
                               sum(data$Sum_m_k[non_zero_theta] * log(theta[non_zero_theta])) - 
                               sum(data$m_t[time_non_zero] * log(normalized_const[time_non_zero])) + 
                               (shape - 1) * (sum(log(f[non_zero_f_temp]))) - rate * sum(f[non_zero_f_temp]) - 
                                sum(lambda * w_k[non_zero_theta[-c(1,num_nonzero)]]*(log(theta[non_zero_theta[-c(1,2)]]) + 
                                   log(theta[non_zero_theta[-c(num_nonzero,num_nonzero - 1)]]) - 
                                    2 * log(theta[non_zero_theta[-c(1,num_nonzero)]]))^2))
        }
        else if (TRUE == only_PA) { 
            log_likelihood <-  c(log_likelihood, sum(data$Sum_m_k[non_zero_theta] * log(theta[non_zero_theta])) - 
                                 sum(data$m_t[time_non_zero] * log(normalized_const[time_non_zero])) -  
                                  sum(lambda * w_k[non_zero_theta[-c(1,num_nonzero)]]*(log(theta[non_zero_theta[-c(1,2)]]) + 
                                  log(theta[non_zero_theta[-c(num_nonzero,num_nonzero - 1)]]) - 
                                  2 * log(theta[non_zero_theta[-c(1,num_nonzero)]]))^2))
        }
        else  log_likelihood    <- c(log_likelihood, sum(data$z_j[non_zero_f_temp] * log(f[non_zero_f_temp])) - 
                                     sum(data$m_t[time_non_zero] * log(normalized_const[time_non_zero])) + 
                                     (shape - 1) * (sum(log(f[non_zero_f_temp]))) - rate * sum(f[non_zero_f_temp]))
       
        ######################### quasi-Newton acceleration #########################
        ####  calculate the smallest approximation to the second derivative at current point ####
        if (q > 1) {
            flag <- 0  
        if (i > q) {
            current_pos <- c(theta[non_zero_theta],f[non_zero_f])
            U <- parameter_save[,1:(q-1)] - parameter_save[,q]
            V <- parameter_save[,2:q]     - current_pos            
                                
            if ((FALSE == only_f) && (FALSE == only_PA)) {
                  candidate <- tryCatch(current_pos + step_size * ifelse(q > 2, V%*%solve(crossprod(U,U) - 
                                       crossprod(U,V),crossprod(U, V[,q-1])), 1/(sum(U* U) - sum(U*V)) * sum(V*U) * V), 
                                       error = function(e) {#print("Problem with inversion"); 
                                                            flag <- 1; return(current_pos)})
                 positive <- prod(candidate > 0)      
            }
            else
            if (TRUE == only_PA) {
                candidate[1:num_nonzero] <- tryCatch(current_pos[1:num_nonzero] + step_size * ifelse(q > 2, V[1:num_nonzero,]%*%solve(crossprod(U[1:num_nonzero,],U[1:num_nonzero,]) - 
                                            crossprod(U[1:num_nonzero,],V[1:num_nonzero,]),crossprod(U[1:num_nonzero,] , V[1:num_nonzero,q-1])), 
                                            1/(sum(U[1:num_nonzero]* U[1:num_nonzero]) - sum(U[1:num_nonzero]*V[1:num_nonzero])) * 
                                              sum(V[1:num_nonzero]*U[1:num_nonzero]) * V[1:num_nonzero]), 
                                            error = function(e) {#print("Problem with inversion"); 
                                                                 flag <- 1 ;return(current_pos[1:num_nonzero])})
                positive <- prod(candidate[1:num_nonzero] > 0) 
            }
            else {
                candidate[-(1:num_nonzero)] <- tryCatch(current_pos[-(1:num_nonzero)]  + 
                                               step_size * ifelse(q > 2, V[-(1:num_nonzero),]%*%solve(crossprod(U[-(1:num_nonzero),],U[-(1:num_nonzero),]) - 
                                               crossprod(U[-(1:num_nonzero),],V[-(1:num_nonzero),]),crossprod(U[-(1:num_nonzero),] ,V[-(1:num_nonzero),q-1])), 1/(sum(U[-(1:num_nonzero)]* U[-(1:num_nonzero)]) - 
                                               sum(U[-(1:num_nonzero)]*V[-(1:num_nonzero)])) * 
                                               sum(V[-(1:num_nonzero)]*U[-(1:num_nonzero)]) * V[-(1:num_nonzero)]),
                                               error = function(e){#print("problem with inversion"); 
                                                                   flag <- 1; return(current_pos[-(1:num_nonzero)])})
                positive <- prod(candidate[-(1:num_nonzero)] > 0) 
            }
            if ((0 == flag) && (1 == positive)) {
                candidate_ok <- candidate_ok + 1
                
                theta_temp                 <- theta
                theta_temp[non_zero_theta] <- candidate[1:num_nonzero]
                f_temp                     <- f
                f_temp[non_zero_f]         <- candidate[-(1:num_nonzero)]
                normalized_const_temp      <- rep(0, dim(data$n_tk)[1])
                if (FALSE == only_PA) {
                    .normalized_constant(normalized_const_temp,data$node_degree,theta_temp,f_temp,data$offset_tk) 
                    
                }  
                else if (TRUE == only_PA)
                    normalized_const_temp <- as.vector(data$n_tk[,non_zero_theta]%*%candidate[1:num_nonzero])
                
                time_temp <- which(normalized_const_temp != 0)   
                if ((FALSE == only_PA) && (FALSE == only_f)) {
                    log_candidate   <- sum(data$z_j[non_zero_f] * log(candidate[-(1:num_nonzero)])) + 
                                           sum(data$Sum_m_k[non_zero_theta] * log(candidate[1:num_nonzero])) - 
                                           sum(data$m_t[time_temp] * log(normalized_const_temp[time_temp])) + 
                                           (shape - 1) * (sum(log(candidate[-(1:num_nonzero)]))) - 
                                           rate * sum(candidate[-(1:num_nonzero)]) - 
                                            sum(lambda * w_k[non_zero_theta[-c(1,num_nonzero)]]*(log(candidate[(1:num_nonzero)[-c(1,2)]]) + 
                                            log(candidate[(1:num_nonzero)[-c(num_nonzero,num_nonzero - 1)]]) - 
                                            2 * log(candidate[(1:num_nonzero)[-c(1,num_nonzero)]]))^2)
                }
                else if (TRUE == only_PA) {
                    log_candidate    <-  sum(data$Sum_m_k[non_zero_theta] * log(candidate[1:num_nonzero])) - 
                                         sum(lambda * w_k[non_zero_theta[-c(1,num_nonzero)]]*(log(candidate[(1:num_nonzero)[-c(1,2)]]) + 
                                         log(candidate[(1:num_nonzero)[-c(num_nonzero,num_nonzero - 1)]]) - 
                                         2 * log(candidate[(1:num_nonzero)[-c(1,num_nonzero)]]))^2)
                }
                else log_candidate   <- sum(data$z_j[non_zero_f] * log(candidate[-(1:num_nonzero)])) - 
                                        sum(data$m_t[time_temp] * log(normalized_const_temp[time_temp])) + 
                                        (shape - 1) * (sum(log(candidate[-(1:num_nonzero)]))) - rate * sum(candidate[-(1:num_nonzero)]) 
                  
               
                
                if (log_candidate > log_likelihood[length(log_likelihood)]) {
                    
                    theta[non_zero_theta] <- candidate[1:num_nonzero]
                    f[non_zero_f]         <- candidate[-(1:num_nonzero)]
                    log_likelihood[length(log_likelihood)]  <- log_candidate
                    normalized_const                        <- normalized_const_temp
                }
                else {
                   #print(paste("ll of rejected candidate:",log_candidate))
                }
            }
            parameter_save[,1:(q-1)]           <- parameter_save[,2:q] 
            parameter_save[1:num_nonzero,q]    <- theta[non_zero_theta]
            parameter_save[-(1:num_nonzero),q] <- f[non_zero_f]  
        }
        else {
            parameter_save[1:num_nonzero,i]    <- theta[non_zero_theta]
            parameter_save[-(1:num_nonzero),i] <- f[non_zero_f]
 
        }
        }
        ############### End of quasi-Newton acceleration ################
        
        if (TRUE == debug)
            print(log_likelihood[length(log_likelihood)])
        break_flag <- FALSE
        if (TRUE == auto_stop)
            if (length(log_likelihood) > 1)
             tryCatch({if (abs(log_likelihood[length(log_likelihood)] - log_likelihood[length(log_likelihood) - 1]) / 
                 (abs(log_likelihood[length(log_likelihood) - 1]) + 1) < stop_cond)
                  break_flag <- TRUE;},error = function(e) { #print(as.vector(normalized_const));print(f[non_zero_f]);
                                                             #print(non_zero_f);
                                                             break_flag <- TRUE;})
        if (break_flag)
            break;   
        ##################### Update f ######################
        if (FALSE == only_PA){
            .update_f(f,non_zero_f,data$node_degree,theta,data$z_j,normalized_const,data$m_t,shape,rate)
            .normalized_constant(normalized_const, data$node_degree,theta,f,data$offset_tk)
            time_non_zero     <- which(normalized_const != 0)
        }
        
        #####################  Update A #######################
        if (FALSE == only_f) {
            if (FALSE == only_PA) {
                #print("problem in coeff theta")
                temp5  <- .coeff_theta(data$node_degree, f, normalized_const,data$m_t,data$start_deg + data$G)     
                temp4  <- temp5[non_zero_theta] + colSums(data$m_t[time_non_zero] / normalized_const[time_non_zero] * 
                                              data$offset[time_non_zero,non_zero_theta])
            }
            else 
                temp4 <- colSums(data$m_t[time_non_zero]/normalized_const[time_non_zero] * 
                                 data$n_tk[time_non_zero,non_zero_theta])
            if (lambda <= 0) {
                theta[non_zero_theta] <- data$Sum_m_k[non_zero_theta]/temp4
            }
            else {
                g_1  <- function(x) {
                        (data$Sum_m_k[non_zero_theta[1]] - 2*w_k[non_zero_theta[2]] * lambda * (-2 * log(theta[non_zero_theta[2]]) + 
                        log(theta[non_zero_theta[3]]) - 3 * log(theta[non_zero_theta[1]])))/x -   
                         temp4[1] - 8 * lambda * w_k[non_zero_theta[2]] * log(x) / x}
                g_2 <- function(x) {
                       (data$Sum_m_k[non_zero_theta[2]] - 2*lambda * 
                          (2*w_k[non_zero_theta[2]]*(-2*log(theta[non_zero_theta[2]]) -log(theta[non_zero_theta[3]]) - log(theta[non_zero_theta[1]])) + 
                             w_k[non_zero_theta[3]]*(-2*log(theta[non_zero_theta[3]]) +log(theta[non_zero_theta[4]]) -3*log(theta[non_zero_theta[2]]) ))) / x - 
                         temp4[2] - (16*w_k[non_zero_theta[2]] + 8*w_k[non_zero_theta[3]]) * lambda * log(x) / x}
               g_semiend <- function(x) {
                            (data$Sum_m_k[non_zero_theta[num_nonzero - 1]] - 2*lambda * 
                            (2*w_k[non_zero_theta[num_nonzero - 1]]*(-2*log(theta[non_zero_theta[num_nonzero - 1]]) - log(theta[non_zero_theta[num_nonzero]]) - log(theta[non_zero_theta[num_nonzero - 2]])) +
                             w_k[non_zero_theta[num_nonzero -2]] *(-3*log(theta[non_zero_theta[num_nonzero - 1]]) + log(theta[non_zero_theta[num_nonzero - 3]]) - 2*log(theta[non_zero_theta[num_nonzero - 2]])))) / x - 
                              temp4[num_nonzero - 1] - 
                              (16*w_k[non_zero_theta[num_nonzero - 1]] + 8*w_k[non_zero_theta[num_nonzero -2]]) * lambda * log(x) / x}
               g_end     <-  function(x) {
                             (data$Sum_m_k[non_zero_theta[num_nonzero]] - 2*lambda * w_k[non_zero_theta[num_nonzero -1]] *(-3*log(theta[non_zero_theta[num_nonzero]]) + log(theta[non_zero_theta[num_nonzero - 2]]) - 2*log(theta[non_zero_theta[num_nonzero - 1]]))) / x  -   
                             temp4[num_nonzero] - 8*w_k[non_zero_theta[num_nonzero -1]]*lambda*log(x)/x}  

               U_1 <- data$Sum_m_k[update_theta] - 2 * lambda * (2*w_k[update_theta]*(-2*log(theta[update_theta]) - log(theta[plus_1]) - log(theta[minus_1])) + 
                                                                   w_k[minus_1]*(-3*log(theta[update_theta]) + log(theta[minus_2]) - 2*log(theta[minus_1])) + 
                                                                   w_k[plus_1]*(-2*log(theta[plus_1]) + log(theta[plus_2]) - 3*log(theta[update_theta])))
               U_2 <- (16*w_k[update_theta] + 8*w_k[minus_1] + 8*w_k[plus_1]) * lambda
               U_3 <- temp4[-c(1,2,num_nonzero-1,num_nonzero)]
               theta[non_zero_theta[1]] <- tryCatch(uniroot(g_1,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                    error = function(e) {return(theta[non_zero_theta[1]])})
        
               theta[non_zero_theta[2]] <- tryCatch(uniroot(g_2,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                   error = function(e) {return(theta[non_zero_theta[2]])})
               #parallelization here
               for (jj in 1:length(update_theta)) {
                  g <- function(x){U_1[jj]/x - U_2[jj]*log(x)/x - U_3[jj]}
                  theta[update_theta[jj]] <- tryCatch(uniroot(g,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                      error = function(e) theta[update_theta[jj]])
               }
               theta[non_zero_theta[num_nonzero - 1]] <- tryCatch(uniroot(g_semiend,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                                  error = function(e) 
                                                                  {return(theta[non_zero_theta[num_nonzero - 1]])})
               theta[non_zero_theta[num_nonzero]] <- tryCatch(uniroot(g_end,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                              error = function(e) {return(theta[non_zero_theta[num_nonzero]])})
             }
        }
}
    ########## End of Iteration ########################    
        
      if (normalized_f == TRUE) {
          f <- length(f) * f / sum(f)
          if (FALSE == only_PA){
            .normalized_constant(normalized_const,data$node_degree,theta,f,data$offset_tk) 
          }
          else normalized_const <- as.vector(data$n_tk[,non_zero_theta]%*%theta[non_zero_theta])
      }
       ##### Variance of f ####################
       cov_f   <- rep(0,length(f))
       if (FALSE == only_PA)
           .cal_var_f(cov_f,non_zero_f,data$node_degree,theta,f,data$z_j,
                     normalized_const,data$m_t,shape)
    
      hessian_of_regularization <- function(theta){
          n      <- length(theta)
          result <- vector()
          if (n - 2 >= 3){
              result <- c(result,2*w_k[non_zero_theta[2]]*((1-log(theta[1]))/theta[1]^2 + (2*log(theta[2]) - log(theta[3])) / theta[1]^2 ) )
              
              result <- c(result,2*w_k[non_zero_theta[3]]*((1-log(theta[2]))/theta[2]^2 + (2*log(theta[3]) - log(theta[4])) / theta[2]^2 ) +
                                 2*w_k[non_zero_theta[2]]*((2-2*log(theta[2]))/theta[2]^2 + (log(theta[3]) + log(theta[1])) / theta[2]^2 ))
              for (ii in 3:(n-2))
                  result <- c(result,2*w_k[non_zero_theta[ii]]*((2-2*log(theta[ii]))/theta[ii]^2 + (log(theta[ii+1]) + log(theta[ii-1])) / theta[ii]^2 ) + 
                                     2*w_k[non_zero_theta[ii + 1]]*((1-log(theta[ii]))/theta[ii]^2 + (2*log(theta[ii+1]) - log(theta[ii+2])) / theta[ii]^2 ) +
                                     2*w_k[non_zero_theta[ii-1]]*((1-log(theta[ii]))/theta[ii]^2 + (2*log(theta[ii-1]) - log(theta[ii-2])) / theta[ii]^2 ))
              ii <- n - 1
              result <- c(result, 2*w_k[non_zero_theta[ii]]*((2-2*log(theta[ii]))/theta[ii]^2 + (log(theta[ii+1]) + log(theta[ii-1])) / theta[ii]^2 ) +
                                  2*w_k[non_zero_theta[ii-1]]*((1-log(theta[ii]))/theta[ii]^2 + (2*log(theta[ii-1]) - log(theta[ii-2])) / theta[ii]^2 ))
              ii <- n
              result <- c(result, 2*w_k[non_zero_theta[ii-1]]*((1-log(theta[ii]))/theta[ii]^2 + (2*log(theta[ii-1]) - log(theta[ii-2])) / theta[ii]^2 )) 
          }
              return(w_k[non_zero_theta]*result)
      }

      #interpolation
    if (TRUE == interpolate) {
    theta_nonzero <- which(theta != 0)
    if (length(theta_nonzero) > 0)
      if (theta_nonzero[1] > 1) {
          theta[1:(theta_nonzero[1] - 1)]   <- theta[theta_nonzero[1]]  
          cov_bin[1:(theta_nonzero[1] - 1)] <- cov_bin[theta_nonzero[1]] 
      }
      if (length(theta_nonzero) > 1) {
          for (i in 1:(length(theta_nonzero) - 1))
              if (theta_nonzero[i+1] > theta_nonzero[i] + 1) {
              regress <- lm(c(log(theta[theta_nonzero[i]]),log(theta[theta_nonzero[i+1]]))~ 
                         c(log(center_k[theta_nonzero[i]]),log(center_k[theta_nonzero[i+1]]))) 
              for (j in (theta_nonzero[i] + 1):(theta_nonzero[i+1] - 1))
                  theta[j] <- exp(log(center_k[j])*regress$coefficients[2] + 
                                    regress$coefficients[1])           
              }  
      }
      if (length(theta_nonzero) > 0)
         if (theta_nonzero[length(theta_nonzero)] < length(theta))
             theta[(theta_nonzero[length(theta_nonzero)] + 1):length(theta)] <- 
              theta[theta_nonzero[length(theta_nonzero)]]  
    }    
    # get the center of each bin
    center_k    <- rep(0, length(theta)) 
    center2_k   <- center_k
    if (data$start_deg > 0) {
        center_k[1:data$start_deg]  <- 0:(data$start_deg - 1)
        center2_k[1:data$start_deg] <- 0:(data$start_deg - 1)
    }
    for (i in 1:data$G) {
        if (data$begin_deg[i] != 0) {
    #              center_k[i]  <- round((data$begin_deg[i] + data$end_deg[i])/2)  
            center_k[data$start_deg + i] <- round(data$begin_deg[i]*sqrt((data$begin_deg[i] + data$interval_length[i] - 1)/ data$begin_deg[i]))
            #center2_k[data$start_deg + i] <- round((data$begin_deg[i] + data$end_deg[i])/2)  
            center2_k[data$start_deg + i] <- data$end_deg[i]
        } else {
              center_k[data$start_deg + i] <- data$end_deg[i]  
              #center2_k[i]  <- round((data$begin_deg[i] + data$end_deg[i])/2)   
              center2_k[data$start_deg + i] <- data$end_deg[i]
          }
  #        
    }
    if (only_PA == FALSE) {
        beg     <- which(center_k >= data$deg_thresh & theta != 0)[1]
    } else beg <- which(theta != 0)[1]    
    theta <- theta/theta[beg]
   

    if (FALSE == only_PA){
        .normalized_constant(normalized_const,data$node_degree,theta,f,data$offset_tk) 
    }
        else normalized_const <- as.vector(data$n_tk[,non_zero_theta]%*%theta[non_zero_theta])
    ############# Calculate the variance of theta #####################
    cov_bin1       <- rep(0,length(theta)) 
    non_zero       <- which(theta != 0)
    non_zero_theta <- non_zero
    time_non_zero     <- which(normalized_const != 0)
    if (FALSE == only_PA) {
      temp4  <- .coeff_var(data$node_degree, f, normalized_const,data$m_t,data$offset_tk, data$start_deg + data$G) 
      temp4 <- temp4[non_zero]
    }
   
    else {
      temp4 <- colSums(data$n_tk[time_non_zero,non_zero_theta]^2 * data$m_t[time_non_zero] / normalized_const[time_non_zero]^2)
    }
    aa <- 1 / (data$Sum_m_k[non_zero_theta] /theta[non_zero_theta] ^ 2 + - temp4 + 
                 lambda * hessian_of_regularization(theta[non_zero_theta]))
    bb <-  1 / (data$Sum_m_k[non_zero_theta] /theta[non_zero_theta] ^ 2 +  
                  lambda * hessian_of_regularization(theta[non_zero_theta])) 
    cov_bin1[non_zero] <- ifelse(aa > 10^-20, aa, bb) 
   
    cov_bin1[cov_bin1 == Inf] <- 0

     cov_bin <- cov_bin1
     var_log_bin  <- cov_bin/ifelse(theta != 0, theta^2,1)
     upper_bin    <- exp(log(theta) + 2*sqrt(var_log_bin))
     lower_bin    <- exp(log(theta) - 2*sqrt(var_log_bin))
     

     non_zero_center <- center2_k > 0 & theta > 0 &  var_log_bin > 0 
     if (FALSE == only_f)
         alpha_center <- lm(log(theta[non_zero_center]) ~ log(center2_k[non_zero_center]),weights = 1/var_log_bin[non_zero_center])$coefficients[2]
     else 
         alpha_center <- NULL   
     ############# Return theta to A #####################################
      A                           <- rep(0,data$deg.max)
      cov                         <- rep(0,data$deg.max)   
      weight_A                    <- rep(1,data$deg.max)
      A[1:(data$start_deg + 1)]        <- theta[1:(data$start_deg + 1)] 
      cov[1:(data$start_deg + 1)]      <- cov_bin[1:(data$start_deg + 1)]
      weight_A[1:(data$start_deg + 1)] <- 1
      for (i in 1:data$G) {
          weight_A[(data$begin_deg[i]:data$end_deg[i]) + 1] <- data$interval_length[i]             
          A[(data$begin_deg[i]:data$end_deg[i]) + 1]        <- theta[data$start_deg + i]
          cov[(data$begin_deg[i]:data$end_deg[i])+1]        <- cov_bin[data$start_deg + i] #* weight_A[data$begin_deg[i] + 1]
      }
      interval     <- 0:(length(A) - 1)
      
      #return(list(k = interval, A = A))
      non_zero     <- which(A > 10^-20 & cov > 10^-20)
      non_zero     <- non_zero[non_zero >= data$deg_thresh]
      k_non_zero   <- interval[non_zero]
    
      A            <- A[non_zero] 
      #cc           <- exp(mean(log(k_non_zero[k_non_zero > 0])) - mean(log(A[A > 0])))
      cc           <- 1
      A            <- cc * A
      weight_A     <- weight_A[non_zero]
      cov          <- cc ^ 2 * cov[non_zero]  
      ############### fitting A_k = k^alpha ##################
      var_log     <- cov / (A ^ 2)
      sd_log      <- sqrt(var_log)
      log_A       <- log(A)
      log_k       <- log(k_non_zero)
      
      if ((only_f == FALSE) && (length(k_non_zero) > 0))  {
              if (0 == k_non_zero[1]) 
                  linear_fit  <- lm(log_A[-1] ~ log_k[-1] , weights = 1 / (weight_A[-1]*var_log[-1])) 
              else
                  linear_fit  <- lm(log_A ~log_k , weights = 1 / (weight_A * var_log)) 
      }
      else linear_fit <- list(coefficients=c(-1,-1))
      names(linear_fit$coefficients) <- c("offset","attachment exponent")
      upper_A       <- exp(log(A) + 2 * sd_log)
      lower_A       <- exp(log(A) - 2 * sd_log)
      
      ##### Variance of f ####################
    
      if (rate != 0)
          f[-non_zero_f] <- shape / rate
      else f[-non_zero_f] <- 0

      f_new        <- rep(1,data$N)
      names(f_new) <- data$node_id

      f_new[as.character(data$f_position)] <- f
      cov_f_new                  <- rep(0,data$N)
      names(cov_f_new) <- data$node_id
      cov_f_new[as.character(data$f_position)] <- cov_f
      non_zero_f          <- f_new > 10^-20 & cov_f_new > 10^-20
    
      upper_f             <- rep(0,data$N)
      upper_f[non_zero_f] <- exp(log(f_new[non_zero_f]) + 2 * sqrt(cov_f_new[non_zero_f] / f_new[non_zero_f] ^ 2))

      lower_f             <- rep(0,data$N)
      lower_f[non_zero_f] <- exp(log(f_new[non_zero_f]) - 2 * sqrt(cov_f_new[non_zero_f] / f_new[non_zero_f] ^ 2))

      result <- list(A           = A   ,          var_A        = cov,      linear_fit     = linear_fit, mode_f = mode_f[1], center_k = center_k,
                     theta = theta, upper_bin = upper_bin, lower_bin = lower_bin,
                   k             = k_non_zero ,   weight_of_A  = weight_A, var_logA       = var_log,  #alpha_center = alpha_center,      
                   upper_A       = upper_A,       lower_A      = lower_A,  #alpha = alpha_center,
                   alpha          = linear_fit$coefficients[2], 
                   var_bin = cov_bin,
                   f             = f_new,         var_f        = cov_f_new, only_PA = only_PA, only_f = only_f,
                   upper_f       = upper_f,       lower_f      = lower_f,  log_likelihood = log_likelihood, lambda = lambda, 
                   shape = shape, rate = rate, normalized_f = normalized_f, deg_threshold = data$deg_thresh,
                   stop_cond = stop_cond, auto_lambda = auto_lambda, ratio = ratio, G = data$G,shape = shape, rate = rate)
      class(result) <- "PAFit"
      return(result)
}
