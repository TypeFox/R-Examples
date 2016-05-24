#~~~~~~ Forecast_comb_all   ~~~~~~~~~~~~~~~~~~~~~~~~~~

Forecast_comb_all <-  function(obs, fhat, fhat_new = NULL){
  TT <-  NROW(fhat)
  l <-  NCOL(fhat)
  pb <- txtProgressBar(min = 1, max= l, style = 3)
  h0 = list()
  MAL0 <- HQ0 <- BIC0 <- AICc0 <- AIC0 <- pred <- list()
  lm_full <- lm(obs ~ as.matrix(fhat))
  sig_full <- summary(lm_full)$sig^2
  aic_full <- as.numeric(-2*logLik(lm_full) + 2*(l+1) )
  bic_full <- as.numeric(-2*logLik(lm_full) + log(TT)*(l+1) )
  aicc_full <- as.numeric(aic_full+ 2*(l+2)*(l+3)/(TT-2) )
  hq_full <- as.numeric(-2*logLik(lm_full) + log(log(TT))*(l+1))
  full_model_crit <- t(as.matrix(c(aic_full, aicc_full, bic_full, hq_full)))
  colnames(full_model_crit) <- c("AIC", "AICc", "BIC", "HQ")
  for (k in 1:l){
    h0[[k]] = as.matrix( combn(l, k) )
    crit <- lm0 <-  list()
    pred0 = matrix(nrow = TT, ncol = NCOL(h0[[k]]))
    for ( j in 1:NCOL(h0[[k]]) ) {
      lm0[[j]] = lm(obs ~ as.matrix(fhat[, h0[[k]][,j] ]))
      crit$AIC[[j]] <- as.numeric(-2*logLik(lm0[[j]]) + 2*(k+1) )
      crit$AICc[[j]] <- crit$AIC[[j]] + 2*(k+2)*(k+3)/(TT-k-1) # k does not include the intercept
      crit$BIC[[j]] <-  as.numeric(-2*logLik(lm0[[j]]) + log(TT)*(k+1)) # BIC
      crit$HQ[[j]] <- as.numeric(-2*logLik(lm0[[j]]) + log(log(TT))*(k+1)) # Hannan Quinn
      temp <- summary(lm0[[j]])
      crit$MAL[[j]] <- temp$sig^2*(2*(k+1)) + sig_full/TT # Mallow
      if(is.null(fhat_new)){
        pred0[,j] = t(t(lm0[[j]]$coef) %*% t(as.matrix(cbind(rep(1, TT), fhat[, h0[[k]][,j]]))))
      } else {
        if (j == 1) { pred0 = matrix(nrow = NROW(fhat_new) , ncol = NCOL(h0[[k]])) }
        pred0[, j] = t(lm0[[j]]$coef %*% t(as.matrix(cbind(rep(1, NROW(fhat_new)), fhat_new[, h0[[k]][,j]])))   )
      } }
    pred[[k]] <- pred0
    AIC0[[k]] <- crit$AIC
    AICc0[[k]] <- crit$AICc
    BIC0[[k]] <- crit$BIC
    HQ0[[k]] <- crit$HQ
    MAL0[[k]] <- crit$MAL
    setTxtProgressBar(pb, k)
  }
  pred <- do.call(cbind, pred)
  # get the weights according to the information criteria
  AIC_weights <- exp(-(1/2)*(unlist(AIC0)- (full_model_crit[1])))
  AIC_weights <- AIC_weights/sum(AIC_weights)
  AICc_weights <- exp(-(1/2)*(unlist(AICc0)- (full_model_crit[2])))
  AICc_weights <- AICc_weights/sum(AICc_weights)
  BIC_weights <- exp(-(1/2)*(unlist(BIC0)- (full_model_crit[3])))
  BIC_weights <- BIC_weights/sum(BIC_weights)
  HQ_weights <- exp(-(1/2)*(unlist(HQ0)- (full_model_crit[4])))
  HQ_weights <- HQ_weights/sum(HQ_weights)
  MAL_weights <- (1/unlist(MAL0))/sum(1/unlist(MAL0))
  return(list(pred = pred, full_model_crit = full_model_crit, aic =AIC_weights,
              aicc = AICc_weights, bic = BIC_weights, hq = HQ_weights, mal = MAL_weights ))
}
