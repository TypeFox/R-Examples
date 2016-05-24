conv.geweke.heidel <-
function (samples_l2_param, X_names, Z_names){
  row_names <- array(NA, dim = c(ncol(samples_l2_param),1))
  for (i in 1:length(X_names))
    for (j in 1:length(Z_names)){
      temp <- (i - 1)*length(Z_names) + j
      row_names[temp] <- paste(X_names[i]," : ", Z_names[j])
    }
  pass_ind <- T
  sol_converge <- geweke.diag(as.mcmc(samples_l2_param),0.2,0.5)
  p_values_converge <- 2*pnorm(-abs(sol_converge$z)) 
  decision_conerge <- array(NA,length(p_values_converge)) 
  decision_conerge[p_values_converge>=0.001] <- 'Passed' 
  decision_conerge[p_values_converge<0.001] <- 'Failed' 
  if (sum(p_values_converge<0.001) > 0) pass_ind <- F
  sol_geweke <- cbind(decision_conerge,round(p_values_converge,4))
  rownames(sol_geweke) <- row_names
  colnames(sol_geweke) <- c('Stationarity test','Convergence p value') 
  
  hddata <- samples_l2_param
  colnames(hddata) <- row_names
  sol_heidel <- heidel.diag(as.mcmc(hddata))  
  sol_heidel[,3] <- round(sol_heidel[,3],4)
  if (sum(sol_heidel[,1]<1) > 0) pass_ind <- F
  sol_heidel[which(sol_heidel[,1]==1),1]<-'Passed'
  sol_heidel[which(sol_heidel[,1]<1),1]<-'Failed'
  
  sol <- list(sol_geweke = sol_geweke, sol_heidel = sol_heidel, pass_ind = pass_ind)
}
