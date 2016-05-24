table.ANCOVA <-
function(samples_l1_param, X, Z){
  if (length(attr(Z, 'varNames')) > 1){
    num_l1_v <- ncol(X)
    num_l2_v <- length(attr(Z, 'varNames')) - 1    
    num_id <- nrow(Z)
    ancova_table <- matrix(NA, nrow = num_l1_v, ncol = num_l2_v+2) 
    rownames(ancova_table) <- colnames(X)
    colnames(ancova_table) <- c(attr(Z, 'varNames')[-1], 'Residuals', 'Total') 
    assign <- attr(Z, 'assign')
    Z_factor_index <- 1:num_l2_v
    n_samples <- nrow(samples_l1_param)
    SS_temp <- 0
    for (i in 1:num_l1_v){
      # y <- array(0, dim = c(num_id, 1))
      #for (n_s in 1:n_samples){
      #  for (j in 1:num_id)
      #    y[j] <- samples_l1_param[n_s,(i-1)*num_id + j]
      #  SS <- ssquares(y, Z, assign, Z_factor_index)
      #  SS_temp <- SS_temp + c(SS$factor_SS, SS$SSE, SS$SS_TO)
      #}
      y <- t(as.matrix(samples_l1_param)[1:n_samples,((i-1)*num_id +1) : (i*num_id)])
      SS <- ssquares(y, Z, assign, Z_factor_index)
      SS_temp <- round(c(SS$factor_SS, SS$SSE, SS$SS_TO), digits = 3)
      SS_res <- array(NA, dim = c(1, length(SS_temp)))
      for (n_sstemp in 1:(length(SS_temp) - 2))
        SS_res[n_sstemp] <- paste(SS_temp[n_sstemp], ' (', round((SS_temp[n_sstemp])/SS_temp[length(SS_temp)]*100, digits = 2), '%)', sep="")
      SS_res[(length(SS_temp)-1):length(SS_temp)] = SS_temp[(length(SS_temp)-1):length(SS_temp)]
      ancova_table[i,] <- SS_res
    }
    return(data.frame(ancova_table))
  }else
    return(NA)
}
