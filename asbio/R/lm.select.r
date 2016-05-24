lm.select <- function(lms, deltaAIC = FALSE){
      models <- matrix(ncol = 1, nrow = length(lms))
      res <- matrix(ncol = 5, nrow = length(lms))
      for(i in 1 : length(lms)){
      res[,1][i] <- AIC(lms[[i]])
          p <- length(coef(lms[[i]])) + 1
          n <- length(fitted(lms[[i]]))
      res[,2][i] <- AIC(lms[[i]]) + (2 * p * (p + 1))/(n - p - 1)
      res[,3][i] <- AIC(lms[[i]], k = log(n))
          MSE <- tail(anova(lms[[1]])$"Mean Sq", 1)
          SSEk <- tail(anova(lms[[i]])$"Sum Sq", 1)
      res[,4][i] <- SSEk/MSE - (n - 2 * p)
      res[,5][i] <- press(lms[[i]])
          tm1 <- attributes(lms[[1]]$terms)$term.labels
          tmi <- attributes(lms[[i]]$terms)$term.labels
          if(deltaAIC == FALSE & any(is.na(match(tmi,tm1)))){
          message("Warning: subsequent models not a subset of first model, Cp values invlaid") }
	models[i] <- as.matrix(format(lms[[i]]$call$formula))
      }
res1 <- data.frame(Model = models, AIC = res[,1], AICc = res[,2], BIC = res[,3], Cp = res[,4], PRESS = res[,5])
if(deltaAIC == TRUE){
deltaAIC <- res[,1] - min(res[,1])
rel.likelihood <- exp(-0.5*deltaAIC) 
res1 <- data.frame(Model = models, deltaAIC = deltaAIC, Rel.likelihood = rel.likelihood, Akaike.weight = rel.likelihood/sum(rel.likelihood)) 
}
res1
}



