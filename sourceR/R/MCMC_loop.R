
doMCMC <- function(data = list(source_data = source_data, human_data = human_data, r = r, prev = prev), num, priors, initials, likelihood_dist, 
                   n_iter, posterior, params_cur, adaptive_matrices, mcmc_params, params_fix, data_names, nc_type) {

  ######################################################################################################
  # Initialise everything ##############################################################################
  ######################################################################################################
  save_i <- seq((mcmc_params$burn_in + 1), n_iter, by = mcmc_params$thin) 
  
  adaptive_a <- TRUE
  print_freq <- max(1, round(n_iter / 10))
  
  source_data <- list() 
  human_data <- list()
  for (t in 1 : num$no_T) {
    source_data[[t]] <- data$source_data[[t]]
    human_data[[t]] <- list()
    for (l in 1 : num$no_L) {
      human_data[[t]][[l]] <- c(data$human_data[[t]][[l]]$human)
    }
  }
  
  r <- data$r
  prev <- data$prev
  
  no_I <- num$no_I
  no_J <- num$no_J
  no_T <- num$no_T
  no_L <- num$no_L
  
  acceptance_matrices <- create_acceptance_rate(num, params_fix, likelihood_dist)
  acceptance = acceptance_matrices$acceptance
  acceptance_rate = acceptance_matrices$acceptance_rate
  
  n_accept = acceptance_matrices$n_accept
  n_update = acceptance_matrices$n_update
  
  tuning_r <- adaptive_matrices$r$tuning_r
  n_update_r <- adaptive_matrices$r$n_update_r
  n_accept_r <- adaptive_matrices$r$n_accept_r 
  
  tuning_a <- adaptive_matrices$a$tuning_a
  n_update_a <- adaptive_matrices$a$n_update_a 
  n_accept_a <- adaptive_matrices$a$n_accept_a
  
  if (likelihood_dist == "nbinom") {
    n_update_d <- adaptive_matrices$d$n_update_d 
    tuning_d <- adaptive_matrices$d$tuning_d
    n_accept_d <- adaptive_matrices$d$n_accept_d
  }

  
  ## calculate the length of the time and location names to pad the output so it lines up
  loclen <- nchar(as.character(data_names$location_ids))
  timlen <- nchar(as.character(data_names$time_ids))
  loc_maxlen <- max(loclen)
  tim_maxlen <- max(timlen)

  
  ## C++ code to calculate lambda i #################################################################### 
  calc_li <- function(no_J, no_I, no_T, no_L, r, a, prev, q) {
    .Call('sourceR_calc_li', PACKAGE = 'sourceR', no_J, no_I, no_T, no_L, r, a, prev, q)
  }
  
  ## C++ code to calculate lambda j ####################################################################
  calc_lj <- function(no_J, no_I, no_T, no_L, r, a, prev, q) {
    .Call('sourceR_calc_lj', PACKAGE = 'sourceR', no_J, no_I, no_T, no_L, r, a, prev, q)
  }
  
  ## C++ code to calculate the log likelihood ##########################################################
  if (likelihood_dist == "pois") {
    like_calc <- function(no_J, no_I, no_T, no_L, lambda_can, d, dat) {
      .Call('sourceR_like_calc_pois', PACKAGE = 'sourceR', no_J, no_I, no_T, no_L, lambda_can, d, dat)
    }
  } else {
    like_calc <- function(no_J, no_I, no_T, no_L, lambda_can, d, dat) {
      .Call('sourceR_like_calc_nbinom', PACKAGE = 'sourceR', no_J, no_I, no_T, no_L, lambda_can, d, dat)
    }
  }
  
  update_ll <- function(params, data = human_data, prev) {
    lambda_can <- calc_li(no_J, no_I, no_T, no_L, r = params$r, a = params$a, prev, q = params$q)
    ll_can <- like_calc(no_J, no_I, no_T, no_L, lambda_can, d = params$d, dat = data)
    return(list(lam_can = lambda_can, ll_can = ll_can))
  }


  #######################################################################################################
  ## Update Functions ###################################################################################
  #######################################################################################################
  
  ## a update ###########################################################################################
  {
    
    ## Propose new a ######################################################################################
    {
      create_propose_a <- function() { 
        # Only update 1 source effect, then rescale the others so that the sum of the source effects is 1
        return( function(t, l, a_can, j) {
              a_can <- params_cur$a
              index <- j
#               print(n_update_a)
#               print(n_accept_a)
              ### adaptive calculations #################################################
              n_update_a[[t]][[l]][index] <<- n_update_a[[t]][[l]][index] + 1  #  we are updating this a value, so add 1 to the number of times it has been updated
              if (isTRUE(all.equal((n_update_a[[t]][[l]][index] %% 50), 0 ))) { # if this parameter has been updated a multiple of 50 times, calculate whether tuning value needs to be changed
                tuning_a[[t]][[l]][index] <<- tuning_a[[t]][[l]][index] + 
                  sign((n_accept_a[[t]][[l]][index] / n_update_a[[t]][[l]][index]) - 0.45) * 
                  min(0.01, n_update_a[[t]][[l]][index] ^ -0.5)
                if (tuning_a[[t]][[l]][index] <= 0) {
                  tuning_a[[t]][[l]][index] <<- 0.00000001 # the tuning value cannot be negative
                }
              }
              
              # print(tuning_a)
              
              ## Proposal Adaptive a code ##################################################################
              pval <- rbern(1, 0.95)
              
              if (isTRUE(all.equal(pval, 1)) & n_update_a[[t]][[l]][index] > 100) {
                a_can[[t]][[l]][index] <- params_cur$a[[t]][[l]][index] + rnorm(1, mean = 0, sd = tuning_a[[t]][[l]][index]) 
                # Use a Normal random walk, rather than log normal random walk because a log normal walk will 
                # propose big positive jumps and small negative jumps. Good for positively skewed distributions.  
                # Ours is not necessarily positvely skewed, just restricted in range.  
                # So use rejection sampling (ie reject any r porposed that is <=0 or >= 1) and a Normal random walk
              } else if (isTRUE(all.equal(pval, 0)) | n_update_a[[t]][[l]][index] <= 100) {
                a_can[[t]][[l]][index] <- params_cur$a[[t]][[l]][index] + rnorm(1, mean = 0, sd = 0.01)
              } 
              
              diff <- (a_can[[t]][[l]][index] - params_cur$a[[t]][[l]][index]) / ((no_J * no_T * no_L) - 1)
              
              for (t1 in 1 : no_T) {
                for (l1 in 1 : no_L) {
                  if ((t1 == t) && (l1 == l)) 
                    a_can[[t1]][[l1]][-index] <- params_cur$a[[t1]][[l1]][-index] - diff # don't change the one we just updated
                  else 
                    a_can[[t1]][[l1]] <- a_can[[t1]][[l1]] - diff
                }
              }

              res_sum_a = sum(sapply(a_can, function(x) sapply(lapply(x, '['), sum, simplify =TRUE), simplify = TRUE))
              a_can[[t]][[l]][index] = a_can[[t]][[l]][index] + (1 - res_sum_a)
              # error was accumulating so the vector no longer summed to 1
              
              # check if any of the a_cans are now smaller than 0 or larger than 1
              auto_reject = any(sapply(a_can, function(x) sapply(lapply(x, '['), function(y) {any(any(y <= 0) | any(y >= 1))}, simplify =TRUE), simplify = TRUE) == TRUE)

              if (auto_reject == TRUE) {
                return(list(auto_reject = TRUE))
              } else {
                return(list(a = a_can, q_ratio = 0, auto_reject = FALSE)) # Normal random walk so no q ratio
              }
            })
      }
      propose_a <- create_propose_a() # returns a and q_ratio containing the proposed a values and their q_ratio for a log normal random walk update
    }
    
    ## Calculate prior a ##################################################################################
    {
      prior_dist_a <- "dirich" #exp
      create_prior_a_calc <- function(prior_dist_a) {
        if (prior_dist_a == "dirich") {
          return(function(prior, a_can, a_cur) {
            prior_cur <- 0 
            prior_can <- 0 
            for (t1 in 1 : no_T) {
              for (l1 in 1 : no_L) {
                prior_cur <- prior_cur + sum((prior - 1) * log(a_cur[[t1]][[l1]]))
                prior_can <- prior_can + sum((prior - 1) * log(a_can[[t1]][[l1]]))
              }
            }
            return(list(prior_cur = prior_cur, prior_can = prior_can))
          })
        } else {
          stop("prior_dist_a must be \"dirich\"")
        }
      }
      
      prior_a_calc <- create_prior_a_calc(prior_dist_a)
    }
    
    create_update_a <- function() {
        return(function(prior, t, l, j) {
          a_can <- params_cur$a

          ####### Propose new a ############################################################################
          props <- propose_a(t, l, a_can, j)
          if (props$auto_reject == FALSE) {
              a_can <- props$a
              priors_val <- prior_a_calc(prior, a_can = a_can, a_cur = params_cur$a)
            
              q_ratio <- props$q_ratio
            ##################################################################################################
            
            ####### Prior a's (Model Part) ###################################################################
            
              prior_cur <- priors_val$prior_cur
              prior_can <- priors_val$prior_can
            ##################################################################################################
            
              update_ll_posterior_can <- update_ll(params = list(a = a_can, r = params_cur$r, d = params_cur$d, q = params_cur$q), data = human_data, prev)
              update_ll_posterior_cur <- update_ll(params = list(a = params_cur$a, r = params_cur$r, d = params_cur$d, q = params_cur$q), data = human_data, prev)
              ll_can <- update_ll_posterior_can$ll_can
            ll_cur <- update_ll_posterior_cur$ll_can
            
            r <- (ll_can + prior_can) - (ll_cur + prior_cur) + q_ratio # Log Normal Proposal Density
            u <- log(runif (1, min = 0, max = 1))
          } else {
            u <- 1 # auto_reject == TRUE, therefore, we must reject
            r <- 0
            #print("Auto-rejecting")
          }
          if (r > u)
          {
            params_cur$a <<- a_can
            params_cur$li <<- update_ll_posterior_can$lam_can
            return(TRUE)
          } else{
            return(FALSE)
          }
        })
    } 
    
    update_a <- create_update_a()
    
  }
  
  ## d update ###########################################################################################
  {
    create_update_d <- function(likelihood_dist) {
      if (likelihood_dist == "pois") {
        return(function(prior, iter_val) {
          return(FALSE)
        })
      } else {
        return(function(prior, iter_val) {
          d_can <- params_cur$d
          
          ### adaptive calculations #################################################
          n_update_d <<- n_update_d + 1 
          if (isTRUE(all.equal((n_update_d %% 50), 0))) {
            tuning_d <<- tuning_d + 
              sign((n_accept_d / n_update_d) - 0.45) * 
              min(0.01, n_update_d ^ -0.5)
            if (tuning_d <= 0) {
              tuning_d <<- 0.00000001
            }
          }
          ##########################################################################
          
          ####### Proposal Adaptive a code ##################################################################
          pval <- rbern(1, 0.95)
          
          if (isTRUE(all.equal(pval, 1)) & iter_val > 100) {#100
            d_can <- exp(rnorm(n = 1, mean = log(params_cur$d), sd = tuning_d)) 
          }else if (isTRUE(all.equal(pval, 0)) | iter_val <= 100) {
            d_can <- exp(rnorm(n = 1, mean = log(params_cur$d), sd = 0.01))
          } 
          
          q_ratio <- log(d_can) - log(params_cur$d)
          
          ####### Prior d's (Model Part) ###################################################################
          prior_cur <- sum(dlnorm(params_cur$d, meanlog = prior[1], sdlog = prior[2], log = TRUE)) 
          prior_can <- sum(dlnorm(d_can, meanlog = prior[1], sdlog = prior[2], log = TRUE)) 
          ##################################################################################################
          
          update_ll_posterior_can <- update_ll(params = list(a = params_cur$a, r = params_cur$r, d = d_can, q = params_cur$q), data = human_data, prev)
          update_ll_posterior_cur <- update_ll(params = list(a = params_cur$a, r = params_cur$r, d = params_cur$d, q = params_cur$q), data = human_data, prev)
          ll_can <- update_ll_posterior_can$ll_can
          ll_cur <- update_ll_posterior_cur$ll_can
          
          post_val <- (ll_can + prior_can) - (ll_cur + prior_cur) + q_ratio
          u <- log(runif (1, min = 0, max = 1))
          
          if (post_val > u) 
          {
            params_cur$d <<- d_can
            n_accept_d <<- n_accept_d + 1
            return(TRUE)
          } else{
            return(FALSE)
          }
        })
      }
    }
    update_d <- create_update_d(likelihood_dist)
  }
  
  ## r update ###########################################################################################
  { 
    create_update_r <- function(r_fix) {
      if (r_fix == FALSE) {
        return(function(iter_val, t, index, prior) {

          r_can <- params_cur$r
          
          ### adaptive calculations #################################################
          n_update_r[[t]][index[1], index[2]] <<- n_update_r[[t]][index[1], index[2]] + 1  #  we are updating this r value, so add 1 to the number of times it has been updated
          if (isTRUE(all.equal((n_update_r[[t]][index[1], index[2]] %% 50), 0))) { # if this parameter has been updated a multiple of 50 times, calculate whether tuning value needs to be changed
            tuning_r[[t]][index[1], index[2]] <<- tuning_r[[t]][index[1], index[2]] + 
              sign((n_accept_r[[t]][index[1], index[2]] / n_update_r[[t]][index[1], index[2]]) - 0.45) * 
              min(0.01, n_update_r[[t]][index[1], index[2]] ^ -0.5)
            if (tuning_r[[t]][index[1], index[2]] <= 0) {
              tuning_r[[t]][index[1], index[2]] <<- 0.00000001 # the tuning value cannot be negative
            }
          }
          
          ## Proposal Adaptive a code ##################################################################
          pval <- rbern(1, 0.95)
          
          if (isTRUE(all.equal(pval, 1)) & iter_val > 100) {
            r_can[[t]][index[1], index[2]] <- params_cur$r[[t]][index[1], index[2]] + rnorm(1, mean = 0, sd = tuning_r[[t]][index[1], index[2]]) 
            # Use a Normal random walk, rather than log normal random walk because a log normal walk will 
            # propose big positive jumps and small negative jumps. Good for positively skewed distributions.  
            # Ours is not necessarily positvely skewed, just restricted in range.  So use rejection sampling (ie reject any r porposed that is <=0 or >= 1) and a Normal random walk
          } else if (isTRUE(all.equal(pval, 0)) | iter_val <= 100) {
            r_can[[t]][index[1], index[2]] <- params_cur$r[[t]][index[1], index[2]] + rnorm(1, mean = 0, sd = 0.01)
          } 
          
          r_can[[t]][-index[1], index[2]] <- params_cur$r[[t]][-index[1], index[2]] - ((r_can[[t]][index[1], index[2]] - params_cur$r[[t]][index[1], index[2]]) / (no_I - 1))
         
          res_sum_r = sum(r_can[[t]][,index[2]])
          # to stop errors accumulating
          r_can[[t]][index[1], index[2]] <- r_can[[t]][index[1], index[2]] + (1 - res_sum_r)
          if (any(r_can[[t]][, index[2]] < 0) | any(r_can[[t]][, index[2]] > 1)) {
            #print("auto-rejecting r")
            return(FALSE)
          }
          # calculates the prior times the likelihood for the source data (multinomial) on the log scale
          prior_m_lik_cur <- sum((source_data[[t]][, index[2]] + prior[[t]][, index[2]] - 1) * log(params_cur$r[[t]][, index[2]]))
          prior_m_lik_can <- sum((source_data[[t]][, index[2]] + prior[[t]][, index[2]] - 1) * log(r_can[[t]][, index[2]]))
          
          q_ratio <- 0 #log(r_can[[t]][index[1], index[2]] / params_cur$r[[t]][index[1], index[2]])
          update_ll_posterior_can <- update_ll(params = list(a = params_cur$a, r = r_can, d = params_cur$d, q = params_cur$q), data = human_data, prev)
          update_ll_posterior_cur <- update_ll(params = list(a = params_cur$a, r = params_cur$r, d = params_cur$d, q = params_cur$q), data = human_data, prev)
          ll_can <- update_ll_posterior_can$ll_can
          ll_cur <- update_ll_posterior_cur$ll_can
          
          r <- (ll_can + prior_m_lik_can) - (ll_cur + prior_m_lik_cur) + q_ratio
          u <- log(runif (1, min = 0, max = 1))
          
          #cat("ll_can: ", ll_can, ", ll_cur: ", ll_cur, ", prior_can: ", prior_m_lik_can, ", prior_cur: ", prior_m_lik_cur, ", q_ratio: ", q_ratio, "\n")

          if (r > u) 
          {
            n_accept_r[[t]][index[1], index[2]] <<- n_accept_r[[t]][index[1], index[2]] + 1
            params_cur$r[[t]] <<- r_can[[t]]
            # params_cur$li <<- update_ll_posterior_can$lam_can
            return(TRUE)
          } else{
            return(FALSE)
          }
        })
      } else {
        return(function(iter_val, t, index, prior) {
          FALSE
        })
      }
    }
    update_r <- create_update_r(params_fix$r)
  }
  
  ## q update ###########################################################################################
  {
    
    params_theta <- function(a, b, k, clusters) {
      yk <- 0
      sums <- 0
      for (t1 in 1 : no_T) {
        for (l1 in 1 : no_L) {
          yk <- yk + sum(human_data[[t1]][[l1]][clusters == k])
          sums <- sums + sum(params_cur$r[[t1]][which(clusters == k),, drop = F] %*% (params_cur$a[[t1]][[l1]] * prev[[t1]]))
        }
      }
      a_star <- a + yk
      b_star <- b + sums
      #             cat("a_star:", a_star, ", b_star:", b_star, ", sums: ", sums, 
      #                 ", yk: ", yk, ", mean = ", a_star / b_star, ", sums*mean: ", sums * a_star / b_star, "\n")
      
      if (a_star <=0 || b_star <= 0) {
        browser()
      } else {
        return(c(a_star, b_star))
      }
    }
    
    sum_y <- function(h, groups, clusters) {
      y_h <- 0
      for (t1 in 1 : no_T) {
        for (l1 in 1 : no_L) {
          if (groups == F) {
            y_h <- y_h + sum(human_data[[t1]][[l1]][h])
          } else {
            y_h <- y_h + sum(human_data[[t1]][[l1]][which(clusters == h)])
          }
          
        }
      }
      return (y_h)
    }
    
    sum_const <- function(h, groups, clusters) {
      sums <- 0
      for (t1 in 1 : no_T) {
        for (l1 in 1 : no_L) {
          if(groups == F) {
            ## calculate const for an individual h
            sums <- sums + sum(params_cur$r[[t1]][h,, drop = F] %*% (params_cur$a[[t1]][[l1]] * prev[[t1]]))
          } else {
            ## calculate sum of consts for all individuals in group h
            sums <- sums + sum(params_cur$r[[t1]][which(clusters == h),, drop = F] %*% (params_cur$a[[t1]][[l1]] * prev[[t1]]))
          }
          
        }
      }
      return (sums)
    }
    
    ########## Update cluster assignments ############
    create_update_cluster <- function(likelihood_dist) {
      if (likelihood_dist == "pois") {
        return(function(prior) {
          
          ## chinese restaurant process
          cluster_can <- params_cur$cluster
          theta_can <- params_cur$theta

          ## for each q, calculate the probability it goes into any of the current groups and the prob it goes in a new group
          for (qi in 1 : no_I) {
            cluster_table <- table(cluster_can)
            n_clusters <- length(cluster_table)
            stopifnot(n_clusters == length(theta_can))
            prob <- numeric(n_clusters + 1)
            
            ## for each q, calculate the probability it goes into any of the current groups 
            ## and the prob it goes in a new group
            ## calculate on log scale and divide by a constant to stop numeric overflow
            
            # log(nk * lambda_i^yi * exp(lambda_i))
            cluster_i <- as.numeric(cluster_can[qi])
            human_i <- sum_y(qi, groups = F, clusters = cluster_can)
            sum_const_i <- sum_const(qi, groups = F, clusters = cluster_can)
            lambda_i <- theta_can * sum_const_i
            
            part_1_a <- log(sum(cluster_table[-cluster_i]))
            part_2_a <- human_i * log(theta_can) - lambda_i
            
            alpha_star <- prior[1] + human_i
            beta_star <- prior[2] + sum_const_i
            
            # log(alpha * (gamma(a_star) / b_star^a_star))
            part_1_b <- log(params_cur$alpha)
            part_2_b <- lgamma(alpha_star) - (alpha_star) * log(beta_star)
            part_3_b <- prior[1] * log(prior[2]) - lgamma(prior[1]) # constant from prior distribution
            
            prob_log <- c(part_1_a + part_2_a, part_1_b + part_2_b + part_3_b) #+ sum_const_i
            
            constant <- max(prob_log)
            prob_log <- prob_log - constant
            prob <- exp(prob_log)

            if(any(!is.finite(prob))) browser()
            
            # choose group for each q wp prob
            new_cluster <- sample(1 : (n_clusters + 1), 1, prob = prob)

            ## if a new group is chosen, draw a theta for it
            if (new_cluster == (n_clusters + 1)) { ## an extra group has been added
              ## add new group to group vec, relevel factor 
              cluster_can <- factor(cluster_can, levels = paste(1 : (n_clusters + 1))) ## add the new factor as a level
              cluster_can[qi] <- new_cluster
              human_i <- sum_y(qi, groups = F)
              sum_const_i <- sum_const(qi, groups = F)
              
              a_star <- prior[1] + human_i # new group, therefore only one y in that group, so don't need to sum
              b_star <- prior[2] + sum_const_i
              new_theta <- Rlab::rgamma(1, a_star, b_star)
              if (new_theta < 0.0000000001) new_theta = 0.0000000001
              theta_can <- c(theta_can, new_theta)
              # if(any(theta <= 0)) browser()
            } else {
              cluster_can[qi] <- new_cluster
            }
            
            # if(any(is.na(cluster_can))) browser()
          
          ## remove any empty groups
          cluster_table <- table(cluster_can)
          empty_clusters <- which(cluster_table ==  0)
          if (length(empty_clusters) != 0) {
            theta_can <- theta_can[-empty_clusters]
            for (t1 in 1 : length(empty_clusters)) {
              clusters_larger <- which(as.numeric(as.character(cluster_can)) > empty_clusters[t1])
              if (length(clusters_larger) != 0) { ## empty groups
                ## the last group is the empty one if no groups larger 
                cluster_can[clusters_larger] <- as.numeric(as.character(cluster_can[clusters_larger])) - 1
              }
            }
            cluster_can <-  factor(cluster_can, levels = paste(1 : length(cluster_table[-empty_clusters]))) 
          }
          if(any(is.na(cluster_can))) browser()
          }
          
          params_cur$theta <<- theta_can
          params_cur$cluster <<- cluster_can
          params_cur$q <<- theta_can[cluster_can]
          
          return (TRUE)
        })
      } else {
        return(function(g_max, prior) {
          FALSE
        })
      }
    }
    
    create_update_theta <- function() {
        return(function(prior) {
          
          theta_can <- params_cur$theta
        
          ### draw new theta's for all clusters
          theta_can <- sapply(1 : length(theta_can), function(k) {
            params <- params_theta(a = prior[1], b = prior[2], k, clusters = params_cur$cluster)
#             print("just before")
#             cat("params[1] ", params[1], ", params[2] ", params[2])
            Rlab::rgamma(1, shape = params[1], rate = params[2])})
          # print("just after")
          
          if (any(theta_can <= 0.0000000001)) {
            theta_can[which(theta_can <= 0.0000000001)] <- 0.0000000001
          }
          
          params_cur$theta <<- theta_can
          params_cur$q <<- theta_can[params_cur$cluster]
          
          return (TRUE)
        })
    }
    
    update_cluster <- create_update_cluster(likelihood_dist)
    update_theta <- create_update_theta()
    
  }

  ## update posterior ###################################################################################
  {
    save_a <- function(params_cur, no_T, no_L) {
      for (t in 1 : no_T) {
        for (l in 1 : no_L) {
          posterior$a[[t]][[l]][k, ] <<- params_cur$a[[t]][[l]]
        }
      }
    }
    
    save_lambdas <- function(params_cur, no_J, no_I, no_T, no_L, k) {
      update_li <- calc_li(no_J, no_I, no_T, no_L, r = params_cur$r, a = params_cur$a, prev = prev, q = params_cur$q) 
      update_lj <- calc_lj(no_J, no_I, no_T, no_L, r = params_cur$r, a = params_cur$a, prev = prev, q = params_cur$q)
      for (t in 1 : no_T) {
        for (l in 1 : no_L) {
          posterior$li[[t]][[l]][k, ] <<- update_li[[t]][[l]]
          posterior$lj[[t]][[l]][k, ] <<- update_lj[[t]][[l]]
        }
      }
    }
    
    
    save_q <- function(params_cur, k) {
      posterior$cluster[k, ] <<- params_cur$cluster
      posterior$theta[k, 1 : length(params_cur$theta)] <<- params_cur$theta
      posterior$q[k, ] <<- params_cur$q
    }
    
    save_q_nbinom <- function(params_cur, k) {
      posterior$q[k, 1] <<- params_cur$q[1]
    }
    
    create_save_r <- function(r_fix) {
      if (r_fix == TRUE) {
        return(function(params_cur, k) {
          FALSE
        })
      } else {
        return(function(params_cur, k) {
          for (t in 1 : no_T) {
            posterior$r[[t]][, , k] <<- params_cur$r[[t]]
          }
        })
      }
    }
    save_r <- create_save_r(params_fix$r)
    
    save_params <- function(params_cur, no_J, no_I, no_T, no_L, k, fix_r, save_lambda, likelihood_dist) {
      save_a(params_cur, no_T, no_L)
      if (save_lambda == TRUE) {
        save_lambdas(params_cur, no_J, no_I, no_T, no_L, k)
      } 
      if (likelihood_dist == "nbinom") {
        posterior$d[k] <<- params_cur$d
        save_q_nbinom(params_cur, k)
      } else {
        save_q(params_cur, k)
      }
      if (fix_r == FALSE) {
        save_r(params_cur, k)
      }
    }
    
  }
  
  ## print acceptance rate
  {
#   print_accept_rate <- function(acceptance, i, no_T, no_L) {
#     if (likelihood_dist == "nbinom") acceptance_rate$d <<- round(acceptance$d / i, digits = 2)
#    
#      for (t in 1 : no_T) {
#        if (params_fix$r == FALSE) acceptance_rate$r[[t]] <<- round(acceptance$r[[t]], digits = 2)
#       for (l in 1 : no_L) {
#         acceptance_rate$a[[t]][[l]] <<- round(acceptance$a[[t]][[l]], digits = 2)
#       }
#     }
#     
#     message("Acceptance rates at iteration ",i,":\n")
#     for (t in 1 : no_T) {
#       if (params_fix$r == FALSE) {
#         message(c("  r              (time ", as.character(data_names$time_ids[t]), ")", paste(rep(" ", tim_maxlen - timlen[t])), 
#                   paste(rep(" ", loc_maxlen - loclen[l])) ,"             :", paste(acceptance_rate$r[[t]], collapse=", ")))
#       }
#       for (l in 1 : no_L) {
#         message(c("  a centered     (time ", as.character(data_names$time_ids[t]), paste(rep(" ", tim_maxlen - timlen[t])), 
#                   ", location ", as.character(data_names$location_ids[l]), paste(rep(" ", loc_maxlen - loclen[l])), ") :", 
#                   paste(acceptance_rate$a[[t]][[l]], collapse=", ")))
#       }
#     }
#     if (likelihood_dist == "nbinom") {
#       message(c("  d: ", acceptance_rate$d))
#     }
#     
#     message("\n")
#   }
  }
  
  ## MCMC Loop ###########################################################################################
  k <- 0

  for (i in 1 : n_iter) {
#     if (isTRUE(all.equal((i %% print_freq), 0))) {
#       if (i > 0) { 
#         print_accept_rate(acceptance, i, no_T, no_L)
#       }
#     }
    
    if (i %in% save_i) {
      k <- k + 1
    }

    ##############################################################################
    ############ Update a ########################################################
    ##############################################################################
    if (prior_dist_a == "dirich") 
      val = no_J
    else 
      val = 1
    for (g in 1 : val) {
      for (t1 in 1 : no_T) { 
        
        for (l1 in 1 : no_L) {
          
          valj <- sample(1 : val, 1)
          
          t <- sample(1 : no_T, 1)
          l <- sample(1: no_L, 1)

          n_update_a[[t]][[l]][valj] <- n_update_a[[t]][[l]][valj] + 1
          n_accept_a[[t]][[l]][valj] <- n_accept_a[[t]][[l]][valj] + 
            update_a(prior = priors$a[[t]][[l]][valj], t = t, l = l, j = valj)
          acceptance$a[[t]][[l]][valj] <- n_accept_a[[t]][[l]][valj] / n_update_a[[t]][[l]][valj]
        }
      }
    }
    
    ##############################################################################
    ############ Update theta and clusters #######################################
    ##############################################################################
    ## all types in 1 group if using negative binomial distribution
    update_cluster(prior = priors$theta) 
    update_theta(prior = priors$theta) 
    
    ##############################################################################
    ############ Update d ########################################################
    ##############################################################################
    acceptance$d <- acceptance$d + update_d(prior = priors$d, iter_val = i)
    
    ##############################################################################
    ############ Update r ########################################################
    ##############################################################################
    for (h in 1 : mcmc_params$n_r) {
      for (t1 in 1 : no_T) {
        t <- sample(1 : no_T, 1)
        update_indices <- c(sample(1 : no_I, 1), sample(1 : no_J, 1))
        n_update$r[[t]] <- n_update$r[[t]] + 1 / mcmc_params$n_r
        n_accept$r[[t]] <- n_accept$r[[t]] + 
          update_r(iter_val = i, t = t, index = update_indices, 
                   prior = priors$r) / mcmc_params$n_r
        acceptance$r[[t]] <- n_accept$r[[t]] / n_update$r[[t]]
      }
    }
    
    ##############################################################################
    ############ Save parameters #################################################
    ##############################################################################
    if (i %in% save_i) {
      save_params(params_cur, no_J, no_I, no_T, no_L, k, params_fix$r, mcmc_params$save_lambda, likelihood_dist)
      save_r(params_cur, k)
    }
  }

#   if (i %% print_freq != 0)
#     print_accept_rate(acceptance, i, no_T, no_L)
  
  source_attrib <- list(posterior = posterior, acceptance = acceptance_rate, data_nested = list(sources = source_data, humans = human_data))
  class(source_attrib) <- "source_attribution"
  return(source_attrib)
}
