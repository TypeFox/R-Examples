library(rstan)
library(mvtnorm)
library(ggplot2)
 
result_fisher_eval = function(nb_params, mat_A_k1, mat_A_k2, n_iter, n_samp, mean_b, params, nb_patients=1){

  Fisher_matrix_covar_temp = array(NA, dim=c(nb_params,nb_params,nb_params,nb_params))
  Fisher_matrix_covar = array(NA, dim=c(nb_params,nb_params,nb_params,nb_params))
  
mean_b = mean_b/(n_samp*2*n_iter)
  Fisher_matrix = crossprod(mat_A_k1[1:n_samp,], mat_A_k2[1:n_samp,])/n_samp
  for(i in 1:nb_params){
    for(j in 1:nb_params){
      for(k in 1:nb_params){
        for(l in 1:nb_params){
          Fisher_matrix_covar_temp[i,j,k,l] = 
          (1/n_samp * sum(mat_A_k1[1:n_samp,i]*mat_A_k2[1:n_samp,j]*mat_A_k1[1:n_samp,k]*mat_A_k2[1:n_samp,l]) - Fisher_matrix[i,j]*Fisher_matrix[k,l])/n_samp
        }
      }
    }
  }
  
  mean_dloglik1 = colMeans(mat_A_k1[1:n_samp,])
  mean_dloglik2 = colMeans(mat_A_k2[1:n_samp,])
  var_dloglik1 = (colMeans(mat_A_k1[1:n_samp,]^2) - mean_dloglik1^2)/n_samp
  var_dloglik2 = (colMeans(mat_A_k2[1:n_samp,]^2) - mean_dloglik2^2)/n_samp
  
  Fisher_matrix = (Fisher_matrix + t(Fisher_matrix))/2
  Fisher_matrix = nb_patients*Fisher_matrix
  if(det(Fisher_matrix) < 1e-02){
    inv_fim = NA
    det_norm_Fisher_matrix = det(Fisher_matrix)^(1/nb_params)
    RSE = NA
    var_det_norm_Fisher_matrix = NA
    print("Warning ! The FIM is singular.")
  }
  else{
    inv_fim = solve(Fisher_matrix)
    det_norm_Fisher_matrix = det(Fisher_matrix)^(1/nb_params)
    RSE = sqrt(diag(inv_fim))/params*100
    var_det_norm_Fisher_matrix = 0
    for(i in 1:nb_params){
      for(j in 1:nb_params){
        for(k in 1:nb_params){
          for(l in 1:nb_params){
            Fisher_matrix_covar[i,j,k,l] = nb_patients^2 *
            (Fisher_matrix_covar_temp[i,j,k,l]+Fisher_matrix_covar_temp[j,i,k,l]+Fisher_matrix_covar_temp[i,j,l,k]+Fisher_matrix_covar_temp[j,i,l,k])/4
            var_det_norm_Fisher_matrix = var_det_norm_Fisher_matrix + 
            (1/nb_params*det_norm_Fisher_matrix)^2 * inv_fim[i,j]*Fisher_matrix_covar[i,j,k,l]*inv_fim[k,l]
          }
        }
      }
    }
  }
  return(list("FIM"=Fisher_matrix,"FIM_covar"=Fisher_matrix_covar, "inv_FIM"=inv_fim, "RSE"=RSE, 
              "det_norm_FIM"=det_norm_Fisher_matrix, "var_det_norm_FIM"=var_det_norm_Fisher_matrix,
               "mean_dloglik1"=mean_dloglik1, "mean_dloglik2"=mean_dloglik2,
               "var_dloglik1"=var_dloglik1, "var_dloglik2"=var_dloglik2, "mean_b"=mean_b))
}


bootstrap_ic_det = function(mat_A_k1, mat_A_k2, n_samp, nb_params, L, normalized=TRUE, nb_patients=1){
  vec_det_fim = numeric(L)
  for(l in 1:L){
    ind_boot = sample(1:n_samp, n_samp, replace=TRUE)
    fim = crossprod(mat_A_k1[ind_boot,], mat_A_k2[ind_boot,])/n_samp
    fim = (fim + t(fim))/2
    fim = nb_patients*fim
    det_fim = det(fim)
    if(normalized==TRUE){
      det_fim = det_fim^(1/nb_params)
    }
    vec_det_fim[l] = det_fim
  }
  binf = quantile(vec_det_fim, 0.025, na.rm=TRUE)
  bsup = quantile(vec_det_fim, 0.975, na.rm=TRUE)
  return(list("binf"=binf, "bsup"=bsup))
}


bootstrap_ic_rse = function(mat_A_k1, mat_A_k2, n_samp, nb_params, L, params, nb_patients=1){
  mat_rse_fim = matrix(0, nrow=L, ncol=nb_params)
  binf = numeric(nb_params)
  bsup = numeric(nb_params)
  for(l in 1:L){
    ind_boot = sample(1:n_samp, n_samp, replace=TRUE)
    fim = crossprod(mat_A_k1[ind_boot,], mat_A_k2[ind_boot,])/n_samp
    fim = (fim + t(fim))/2
    fim = nb_patients*fim
    rse = sqrt(diag(solve(fim)))/params*100
    mat_rse_fim[l,] = rse
  }
  for(p in 1:nb_params){
    binf[p] = quantile(mat_rse_fim[,p], 0.025, na.rm=TRUE)
    bsup[p] = quantile(mat_rse_fim[,p], 0.975, na.rm=TRUE)
  }
  return(list("binf"=binf, "bsup"=bsup))
}
                                            
fisher_evaluation = function(t, y_ini=1, model, model2, model3, params, dim_b, set_seed=TRUE, seed=42, n_samp, n_rep=1, n_iter, n_burn, CV=FALSE, plot_graph=0, L_boot=1000, nb_patients=1){
  if(set_seed==TRUE){
    set.seed(seed)
  }
  nb_t = length(t)
  nb_params = length(params)
mean_b=numeric(dim_b)
  plot_det = c()
  plot_var_det = c()
  plot_binf_boot = c()
  plot_bsup_boot = c()
  if( sum(plot_graph != c(1,2,3,4)) == 4){
    plot_graph = 0
  } 
  
  # Sampling y in its marginal distribution
  sample_y = matrix(NA,nrow=n_samp*n_rep, ncol=nb_t)
  data_cur = list(params=params, mu_b=rep(0, dim_b), t=t, nb_t=nb_t, dim_b=dim_b, n_rep=n_rep)
  if(dim_b > 1){
    init_cur  = list(y=y_ini, b=rmvnorm(1, mean = rep(0, dim_b), sigma = diag(rep(0.1, dim_b)))[1,])
  }
  else if(dim_b == 1){
    init_cur  = list(y=y_ini, b=rnorm(1, 0, sqrt(0.1)))
  }
  temp_sample_y = extract(sampling(model3, data=data_cur, chains=1, init=list(init_cur), warmup=0, iter=n_samp, thin=1, refresh=-1,algorithm="Fixed_param"), permuted=FALSE)
  for(ind_y_samp in 1:n_samp){
    y_samp = temp_sample_y[ind_y_samp,1,which(grepl("y",names(temp_sample_y[1,,]))==TRUE)]
    for(tim in 1:nb_t){
      sample_y[(n_rep*(ind_y_samp-1)+1):(n_rep*ind_y_samp),tim] = y_samp[(n_rep*(tim-1)+1):(n_rep*tim)]
    }
  }

  mat_A_k1 = matrix(NA, nrow=n_samp, ncol=nb_params)
  mat_A_k2 = matrix(NA, nrow=n_samp, ncol=nb_params)
  
  for(ind_y_samp in 1:n_samp){
print(ind_y_samp)
    smp_y = sample_y[(n_rep*(ind_y_samp-1)+1):(n_rep*ind_y_samp),]
    data_cur = list(y=smp_y, params=params, mu_b=rep(0, dim_b), t=t, nb_t=nb_t, dim_b=dim_b, n_rep=n_rep)
    init_cur = list(b=temp_sample_y[ind_y_samp,1,1:dim_b])
    if(dim_b > 1){
      init_cur2 = list(b=rmvnorm(1, mean = rep(0, dim_b), sigma = diag(rep(0.1, dim_b)))[1,])
    }
    else if(dim_b == 1){
      init_cur2 = list(b=rnorm(1, 0, sqrt(0.1)))
    }
    sample_b_sY = extract(sampling(model, data=data_cur, chains=2, init=list(init_cur, init_cur2), warmup=n_burn, iter=10*n_iter+n_burn, thin=10, refresh=-1), permuted=FALSE)
   
    matA_b_k1 = matrix(NA, nrow=n_iter, ncol=nb_params)
    matA_b_k2 = matrix(NA, nrow=n_iter, ncol=nb_params)
    data_cur2 = list(y=smp_y, mu_b=rep(0, dim_b), t=t, nb_t=nb_t, dim_b=dim_b, n_rep=n_rep)
    init_cur2 = list(params=params, b=init_cur$b)
    log_lik = sampling(model2, data=data_cur2, chains=1, init=list(init_cur2), warmup=0, iter=1, refresh=-1, algorithm="Fixed_param")
    for(ind_b in 1:n_iter){
      b_samp1 = sample_b_sY[ind_b,1,1:dim_b]
      b_samp2 = sample_b_sY[ind_b,2,1:dim_b]   
mean_b = mean_b+b_samp1+b_samp2  
      upars1 = unconstrain_pars(log_lik, list(params=params, b=b_samp1))
      upars2 = unconstrain_pars(log_lik, list(params=params, b=b_samp2))
      matA_b_k1[ind_b,] = grad_log_prob(log_lik, upars1, adjust_transform=FALSE)[1:nb_params]
      matA_b_k2[ind_b,] = grad_log_prob(log_lik, upars2, adjust_transform=FALSE)[1:nb_params] 
    }
    mat_A_k1[ind_y_samp,] = colMeans(matA_b_k1[1:n_iter,])
    mat_A_k2[ind_y_samp,] = colMeans(matA_b_k2[1:n_iter,])
    
    if(ind_y_samp>=50 && plot_graph != 0 && ind_y_samp %% 10 == 0){
      res_temp = result_fisher_eval(nb_params, mat_A_k1, mat_A_k2, n_iter, ind_y_samp, mean_b, params, nb_patients)
      plot_det = c(plot_det, res_temp$det_norm_FIM)
      # IC normal
      plot_var_det = c(plot_var_det, res_temp$var_det_norm_FIM)
      plot_interv_inf = pmax(plot_det - 1.96*sqrt(plot_var_det), 0)
      plot_interv_sup = plot_det + 1.96*sqrt(plot_var_det)
      # IC bootstrap
      if(plot_graph == 3 || plot_graph == 4){
        born_boot = bootstrap_ic_det(mat_A_k1, mat_A_k2, ind_y_samp, nb_params, L_boot, normalized=TRUE, nb_patients)
        plot_binf_boot = c(plot_binf_boot, born_boot$binf)
        plot_bsup_boot = c(plot_bsup_boot, born_boot$bsup)      
      }
      lim_y = c(min(plot_interv_inf, plot_binf_boot, na.rm=TRUE), 1.1*max(plot_interv_sup, plot_bsup_boot, na.rm=TRUE))
      x = seq(50, ind_y_samp, by=10) 
      plot(x, plot_det, xlim=c(1,n_samp), xlab="Number of MC samples", ylab="Normalized determinant of the FIM",
      ylim=lim_y, type = "l", col=1, bty='n', lwd=2, main= expression(det(FIM)^frac(1,p)))
      if(plot_graph == 2 || plot_graph == 4){
        lines(x, plot_interv_inf, type = "l", col=3, lty=2)
        lines(x, plot_interv_sup, type = "l", col=3, lty=2)
      }
      if(plot_graph == 3 || plot_graph == 4){
        lines(x, plot_binf_boot, type = "l", col=2, lty=2)
        lines(x, plot_bsup_boot, type = "l", col=2, lty=2)
      }
      if(plot_graph == 2){
        legend(x=0.6*n_samp,y=lim_y[2],legend=c("IC normal"),col=c(3),lty=c(2), cex=1.0,bty='n')
      }
      if(plot_graph == 3){
        legend(x=0.6*n_samp,y=lim_y[2],legend=c("IC bootstrap"),col=c(2),lty=c(2), cex=1.0,bty='n')
      }
      if(plot_graph == 4){
        legend(x=0.6*n_samp,y=lim_y[2],legend=c("IC normal", "IC bootstrap"),col=c(3,2),lty=c(2,2), cex=1.0,bty='n')
      }
    }
  }
  
  res_final = result_fisher_eval(nb_params, mat_A_k1, mat_A_k2, n_iter, n_samp, mean_b, params, nb_patients)
  Fisher_matrix = res_final$FIM
  inv_FIM = res_final$inv_FIM
  Fisher_matrix_covar = res_final$FIM_covar
  det_norm_Fisher_matrix= res_final$det_norm_FIM
  var_det_norm_Fisher_matrix= res_final$var_det_norm_FIM
  born_normal_inf = max(det_norm_Fisher_matrix - 1.96*sqrt(var_det_norm_Fisher_matrix), 0)
  born_normal_sup = det_norm_Fisher_matrix + 1.96*sqrt(var_det_norm_Fisher_matrix)
  born_boot_final = bootstrap_ic_det(mat_A_k1, mat_A_k2, n_samp, nb_params, L_boot, normalized=TRUE, nb_patients)
  mean_dloglik1 = res_final$mean_dloglik1
  mean_dloglik2 = res_final$mean_dloglik2
  var_dloglik1 = res_final$var_dloglik1
  var_dloglik2 = res_final$var_dloglik1
  mean_b = res_final$mean_b
  RSE = res_final$RSE
  boot_rse = bootstrap_ic_rse(mat_A_k1, mat_A_k2, n_samp, nb_params, L_boot, params, nb_patients)
  rse_inf_boot = boot_rse$binf
  rse_sup_boot = boot_rse$bsup
  
  res = NA
  if(CV==TRUE){
    res = list("FIM"=Fisher_matrix, "FIM_covar" = Fisher_matrix_covar, "inv_FIM"=inv_FIM,
               "RSE"=RSE, "RSE_inf_boot"=rse_inf_boot, "RSE_sup_boot"=rse_sup_boot,
               "det_norm_FIM"=det_norm_Fisher_matrix, 
               "det_IC_normal"=c(born_normal_inf, born_normal_sup),
               "det_IC_bootstrap"=c(born_boot_final$binf, born_boot_final$bsup),
               "mean_dloglik1"=mean_dloglik1, "mean_dloglik2"=mean_dloglik2,
               "var_dloglik1"=var_dloglik1, "var_dloglik2"=var_dloglik2,
               "mean_b"=mean_b, "mat_A_k1"=mat_A_k1, "mat_A_k2"=mat_A_k2)
  }
  else{
    res = list("FIM"=Fisher_matrix, "FIM_covar"=Fisher_matrix_covar, "inv_FIM"=inv_FIM,
    "RSE"=RSE, "RSE_inf_boot"=rse_inf_boot, "RSE_sup_boot"=rse_sup_boot, 
    "det_norm_FIM"=det_norm_Fisher_matrix, 
    "det_IC_normal"=c(born_normal_inf, born_normal_sup), 
    "det_IC_boot"=c(born_boot_final$binf, born_boot_final$bsup))
  }
  return(res)
}





template_model = function(path=getwd(), dloglik, nb_t, outcome, nb_params, ind_RE, Cov_list=list(), Sigma_b=FALSE, n_rep=1, name){
  nb_cov = length(Cov_list)
  con = file(paste(path, "/model_", name, ifelse(dloglik==FALSE,1,2), ".stan", sep=""), open = "w")
  dim_b = length(ind_RE)
  
  # data
  cat("data {\n", file = con)
  cat("  int dim_b;\n", file = con)
  if(dim_b > 1){
    cat("  vector[dim_b] mu_b;\n", file = con)
  }
  else{
     cat("  real mu_b;\n", file = con)
  }
  if(outcome=="continuous" || outcome=="binary" || outcome=="longitudinal_binary"|| outcome=="count"){
    cat("  int nb_t;\n", file = con)
    if(nb_t > 1){
      cat("  vector[nb_t] t;\n", file = con)
    }
    else{
      cat("  real t;\n", file = con)
    }
  }
  if(outcome!="continuous" && n_rep > 1){
    cat("  int n_rep;\n", file = con)
  }  
  if(nb_t > 1){
    if(n_rep > 1){
      cat("  int y[n_rep, nb_t];\n", file = con)
    }
    else{
      cat("  vector[nb_t] y;\n", file = con)
    }
  }
  else{
    if(n_rep > 1){
      cat("  vector[n_rep] y;\n", file = con)
    }
    else{
      cat("  real y;\n", file = con)
    }
  } 
  if(dloglik==FALSE){
    cat("  vector[", nb_params,"] params;\n", file = con, sep = "")
  }
  cat("}\n", file = con)
  
  # transformed data
  if(dloglik==FALSE){
    cat("transformed data {\n", file = con)
    if(dim_b > 1){
      if(nb_cov > 0){
        cat("  matrix[dim_b, dim_b] Omega;\n", file = con)
      }
      else{
        cat("  vector[dim_b] Omega;\n", file = con)
      }
    }
    else{
       cat("  real Omega;\n", file = con)
    }
    if(outcome=="continuous" && Sigma_b==FALSE){
      if(nb_t > 1){
        cat("  vector[nb_t] SdEps;\n", file = con)
      }
      else{
        cat("  real SdEps;\n", file = con)
      }
    }
    if(dim_b > 1){
      if(nb_cov > 0){
        cat("  Omega <- diag_matrix(segment(params, ", ind_RE[1], ", dim_b));\n", file = con, sep="")
        for(icov in 1:nb_cov){
          cat("  Omega[", Cov_list[[icov]][1], ", ", Cov_list[[icov]][2], "] <- params[", Cov_list[[icov]][3], "];\n", file = con, sep="")
          cat("  Omega[", Cov_list[[icov]][2], ", ", Cov_list[[icov]][1], "] <- params[", Cov_list[[icov]][3], "];\n", file = con, sep="")
        }
      }
      else{
        cat("  for(i in 1:dim_b){\n", file = con)
        cat("    Omega <- sqrt(params[", ind_RE[1], "+i-1]);\n", file = con, sep="")
        cat("  }\n", file = con)
      }
    }
    else{
       cat("  Omega <- sqrt(params[", ind_RE,"])\n", file = con, sep="")
    }
    if(outcome=="continuous" && Sigma_b==FALSE){
      cat("  SdEps <- TODO;\n", file = con)
    }
    cat("}\n", file = con)
  }
  
  # parameters
  cat("parameters {\n", file = con)
  if(dloglik==TRUE){
    cat("  vector[", nb_params,"] params;\n", file = con, sep = "")
  }
  if(dim_b > 1){
    cat("  vector[dim_b] b;\n", file = con)
  }
  else{
     cat("  real b;\n", file = con)
  }
  cat("}\n", file = con)
  
  # transformed parameters
  if(dloglik==TRUE || Sigma_b==TRUE){
    cat("transformed parameters {\n" ,file = con)
  }
  if(dloglik==TRUE){
    if(dim_b > 1){
      if(nb_cov > 0){
        cat("  matrix[dim_b, dim_b] Omega;\n", file = con)
      }
      else{
        cat("  vector[dim_b] Omega;\n", file = con)
      }
    }
    else{
      cat("  real Omega;\n", file = con)
    }
  }
  if(outcome=="continuous" && (dloglik==TRUE || Sigma_b==TRUE)){
    if(nb_t > 1){
      cat("  vector[nb_t] SdEps;\n", file = con)
    }
    else{
      cat("  real SdEps;\n", file = con)
    }
  }  
  if(dloglik==TRUE){
    if(dim_b >1){
      if(nb_cov > 0){
        cat("  Omega <- diag_matrix(segment(params, ", ind_RE[1], ", dim_b));\n", file = con, sep="") 
        for(icov in 1:nb_cov){
          cat("  Omega[", Cov_list[[icov]][1], ", ", Cov_list[[icov]][2], "] <- params[", Cov_list[[icov]][3], "];\n", file = con, sep="")
          cat("  Omega[", Cov_list[[icov]][2], ", ", Cov_list[[icov]][1], "] <- params[", Cov_list[[icov]][3], "];\n", file = con, sep="")
        }
      }
      else{
        cat("  for(i in 1:dim_b){\n", file = con)
        cat("    Omega <- sqrt(params[", ind_RE[1], "+i-1]);\n", file = con, sep="")
        cat("  }\n", file = con)
      }
    }
    else{
      cat("  Omega <- sqrt(params[", ind_RE, "])\n", file = con, sep="")
    }
  } 
  if(outcome=="continuous" && (dloglik==TRUE || Sigma_b==TRUE)){
    cat("  SdEps <- TODO;\n", file = con)
  }  
  if(dloglik==TRUE || Sigma_b==TRUE){
    cat("}\n" ,file = con)
  }
  
  # model
  cat("model {\n" ,file = con)
  if(nb_cov == 0){
    cat("  b ~ normal(mu_b, Omega);\n", file = con)
  }
  else{
    cat("  b ~ multi_normal(mu_b, Omega);\n", file = con)
  }
  type_link = ifelse(outcome=="continuous", "normal", ifelse(outcome=="binary" || outcome=="longitudinal_binary", "bernoulli_logit", 
              ifelse(outcome=="count", "poisson_log", ifelse(outcome=="time_to_event", "exponential", NA))))
  mod = "  y"
  if(n_rep > 1){
    cat("  for(r in 1:n_rep){\n", file = con)
    mod = paste("  ", mod, "[r]", sep="")
  }
  mod = paste(mod, " ~ ", type_link, "( TODO ", sep="")
  if(outcome=="continuous"){
    mod = paste(mod, ", SdEps);", sep="")
  }
  else{
     mod = paste(mod, ");", sep="")
  }
  cat(mod, "\n", file = con, sep="") 
  if(n_rep > 1){
    cat("  }\n" ,file = con)
  }
  cat("}\n" ,file = con)
  close(con)
}


eval_comb = function(mat_comb, y_ini, model, model2, model3, params, dim_b, set_seed=TRUE, seed, n_samp_mc, n_rep=1, n_iter, n_burn, res_A_prev, L_boot=1000, plot_graph=TRUE, nb_patients=1){
  nb_poss = ncol(mat_comb)
  nb_params = length(params)
  res_det = numeric(nb_poss)
  born_boot_inf = numeric(nb_poss)
  born_boot_sup = numeric(nb_poss)
  fim = vector("list", nb_poss)
  res_A = vector("list", nb_poss)
  for(i_poss in 1:nb_poss){
    res_A[[i_poss]] = vector("list", 2)
    
    eval_fim = fisher_evaluation(t=c(mat_comb[,i_poss]), y_ini=y_ini, model=model, model2=model2, model3=model3, params=params, dim_b=dim_b, set_seed=set_seed, seed=seed, n_samp=n_samp_mc, n_rep=n_rep, n_iter=n_iter, n_burn=n_burn, CV=TRUE, plot_graph = 0, nb_patients=nb_patients)
    if(is.null(res_A_prev)==FALSE){
      mat_A_k1 = rbind(eval_fim$mat_A_k1, res_A_prev[[i_poss]][[1]])
      mat_A_k2 = rbind(eval_fim$mat_A_k2, res_A_prev[[i_poss]][[2]])
    }
    else{
      mat_A_k1 = eval_fim$mat_A_k1
      mat_A_k2 = eval_fim$mat_A_k2
    }
    res_A[[i_poss]][[1]] = mat_A_k1
    res_A[[i_poss]][[2]] = mat_A_k2
    nb_mc = nrow(mat_A_k1)
    fim_temp = crossprod(mat_A_k1, mat_A_k2)/nb_mc
    fim[[i_poss]] = (fim_temp + t(fim_temp))/2 * nb_patients
    res_det[i_poss] = det(fim[[i_poss]])
    bootstrap = bootstrap_ic_det(mat_A_k1, mat_A_k2, nb_mc, nb_params, L_boot, normalized=FALSE, nb_patients)
    born_boot_inf[i_poss] = bootstrap$binf
    born_boot_sup[i_poss] = bootstrap$bsup
  }
  if(plot_graph==TRUE){
    x=1:nb_poss
    y=res_det
    qplot(x,y)+geom_errorbar(aes(x=x, ymin=born_boot_inf, ymax=born_boot_sup), width=0.25)
    ggsave(paste("det_",nb_mc,".pdf",sep=""))
  }
  return(list(res_det, born_boot_sup, res_A, born_boot_inf, fim))
}


phase_select = function(mat_comb, y_ini, model, model2, model3, params, dim_b, set_seed=TRUE, seed, n_samp_mc, n_rep=1, n_iter, n_burn, res_A_prev, L_boot=1000, plot_graph=TRUE, nb_patients=1){
  nb_poss = ncol(mat_comb)
  pre_res = eval_comb(mat_comb=mat_comb, y_ini=y_ini, model=model, model2=model2, model3=model3, params=params,
                       dim_b=dim_b, set_seed=TRUE, seed=seed, n_samp_mc=n_samp_mc, n_rep=n_rep, n_iter=n_iter, 
                       n_burn=n_burn, res_A_prev=res_A_prev, L_boot=L_boot, plot_graph=plot_graph, nb_patients=nb_patients)
  res_det = pre_res[[1]]
  born_boot_sup = pre_res[[2]]
  ind_max = which.max(res_det)
  max_det = res_det[ind_max]
  ind_keep = which(born_boot_sup >= max_det)
  nb_keep = length(ind_keep)
  mat_comb = mat_comb[,ind_keep]
  res_A_prev = vector("list", nb_keep) 
  for(s in 1:nb_keep){
    res_A_prev[[s]] = pre_res[[3]][[ind_keep[s]]]
  }
  return(list(mat_comb, res_A_prev, pre_res))
}


fisher_optimization = function(nb_t, set_t, y_ini, model, model2, model3, params, dim_b, set_seed=TRUE, seed=42, 
step_mc, n_samp_min=30, n_samp_max, n_rep=1, n_iter, n_burn, L_boot=1000, plot_graph=TRUE, nb_patients=1){
  mat_comb = combn(set_t, nb_t)
  nb_poss = ncol(mat_comb)
  nb_params = length(params)
  list_select = vector("list", 1)
  list_select[[1]] = mat_comb
  list_fim = vector("list", 1)
  list_det = vector("list", 1)
  list_binf = vector("list", 1)
  list_bsup = vector("list", 1)
  res_A_prev=NULL
  mc_performed = 0
  p = 0
  res_A = vector("list", 1)
  
  while(ifelse(is.null(ncol(mat_comb)), FALSE, (ncol(mat_comb)>1)) && mc_performed <= n_samp_max){
    p = p+1
    if(mc_performed == 0){
      res = phase_select(mat_comb, y_ini, model, model2, model3, params, dim_b, set_seed=TRUE, seed=seed+mc_performed,
                          n_samp_mc=n_samp_min, n_rep=n_rep, n_iter=n_iter, n_burn=n_burn, 
                          res_A_prev=res_A_prev, L_boot=L_boot, plot_graph=plot_graph, nb_patients=nb_patients)
      mc_performed = mc_performed + n_samp_min 
    }
    else{         
      res = phase_select(mat_comb, y_ini, model, model2, model3, params, dim_b, set_seed=TRUE, seed=seed+mc_performed,
                          n_samp_mc=step_mc, n_rep=n_rep, n_iter=n_iter, n_burn=n_burn, 
                          res_A_prev=res_A_prev, L_boot=L_boot, plot_graph=plot_graph, nb_patients=nb_patients)
      mc_performed = mc_performed + step_mc
    }
    
    res_A_prev = res[[2]]
    mat_comb = res[[1]]
    list_select[[p+1]] = mat_comb
    list_fim[[p]] = res[[3]][[5]]
    list_det[[p]] = res[[3]][[1]]
    list_binf[[p]] = res[[3]][[4]]
    list_bsup[[p]] = res[[3]][[2]]
    res_A[[p]] = res[[3]][[3]]
  }
  if(ncol(mat_comb)>1){
    p = p+1
    ind_max = which.max(res[[3]][[1]])
    mat_comb = mat_comb[,ind_max]
    list_select[[p+1]] = mat_comb
    list_fim[[p]] = res[[3]][[5]][[ind_max]]
    list_det[[p]] = res[[3]][[1]][ind_max]
    list_binf[[p]] = res[[3]][[4]][ind_max]
    list_bsup[[p]] = res[[3]][[2]][ind_max]
    res_A[[p]] = res[[3]][[3]][[ind_max]]
  }
  
  opt_t = mat_comb
  plot_det = c()
  plot_var_det = c()
  plot_binf_boot = c()
  plot_bsup_boot = c()
  for(step_p in 1:(p-1)){
    ind_opt_p = which(list_select[[step_p]][1:nb_t,]==opt_t)
    res_temp = result_fisher_eval(nb_params, res_A[[step_p]][[ind_opt_p]][[1]], res_A[[step_p]][[ind_opt_p]][[2]], n_iter, n_samp_min+step_mc*(step_p-1), 0, params, nb_patients)
    Fisher_matrix = list_fim[[step_p]][[ind_opt_p]]
    Fisher_matrix_covar = res_temp$FIM_covar
    det_norm_Fisher_matrix = list_det[[step_p]][[ind_opt_p]]^(1/nb_params)
    var_det_norm_Fisher_matrix= res_temp$var_det_norm_FIM
    born_normal_inf = det_norm_Fisher_matrix - 1.96*sqrt(var_det_norm_Fisher_matrix)
    born_normal_sup = det_norm_Fisher_matrix + 1.96*sqrt(var_det_norm_Fisher_matrix)
    plot_det = c(plot_det, det_norm_Fisher_matrix)
    # IC normal
    plot_var_det = c(plot_var_det, res_temp$var_det_norm_FIM)
    # IC bootstrap
    born_boot_inf = list_binf[[step_p]][ind_opt_p]
    born_boot_sup = list_bsup[[step_p]][ind_opt_p]
    plot_binf_boot = c(plot_binf_boot, born_boot_inf)
    plot_bsup_boot = c(plot_bsup_boot, born_boot_sup)   
  }
  plot_interv_inf = pmax(plot_det - 1.96*sqrt(plot_var_det), 0)
  plot_interv_sup = plot_det + 1.96*sqrt(plot_var_det)
  lim_y = c(min(plot_interv_inf, plot_binf_boot, na.rm=TRUE), 1.1*max(plot_interv_sup, plot_bsup_boot, na.rm=TRUE))
  plot(n_samp_min+step_mc*(0:(p-2)), plot_det, xlim=c(1,n_samp_min+step_mc(p-2)), xlab="Number of MC samples", ylab="Normalized determinant of the FIM",
  ylim=lim_y, type = "l", col=1, bty='n', lwd=2, main = paste(expression(det(FIM)^frac(1,p)), " for optimal design", sep=""))
  lines(n_samp_min+step_mc*(0:(p-2)), plot_interv_inf, type = "l", col=3, lty=2)
  lines(n_samp_min+step_mc*(0:(p-2)), plot_interv_sup, type = "l", col=3, lty=2)
  lines(n_samp_min+step_mc*(0:(p-2)), plot_binf_boot, type = "l", col=2, lty=2)
  lines(n_samp_min+step_mc*(0:(p-2)), plot_bsup_boot, type = "l", col=2, lty=2)
  legend(x=0.6*n_samp_min+step_mc(p-2),y=lim_y[2],legend=c("IC normal", "IC bootstrap"),col=c(3,2),lty=c(2,2), cex=1.0,bty='n')
  RSE = res_temp$RSE
  boot_rse = bootstrap_ic_rse(res_A[[p-1]][[ind_opt_p]][[1]], res_A[[p-1]][[ind_opt_p]][[2]], n_samp_min+step_mc*(p-2), nb_params, L_boot, params, nb_patients)
  rse_inf_boot = boot_rse$binf
  rse_sup_boot = boot_rse$bsup
  
  res_return = list("opt_t"=mat_comb, "FIM_opt_t"=Fisher_matrix, "FIM_covar_opt_t" = Fisher_matrix_covar,
               "inv_FIM_opt_t" =solve(Fisher_matrix),
               "RSE_opt_t"=RSE, "RSE_inf_boot_opt_t"=rse_inf_boot, "RSE_sup_boot_opt_t"=rse_sup_boot,
               "det_norm_FIM_opt_t"=det_norm_Fisher_matrix, 
               "IC_normal_opt_t"=c(born_normal_inf, born_normal_sup),
               "IC_boot_opt_t"=c(born_boot_inf, born_boot_sup),
               "list_select"=list_select, "list_det"=list_det, 
               "list_boot_inf"=list_binf, "list_boot_sup"=list_bsup, "list_fim"=list_fim)
  
  return(res_return)
}







