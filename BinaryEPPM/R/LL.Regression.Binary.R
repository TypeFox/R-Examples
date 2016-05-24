LL.Regression.Binary <-
function(parameter,model.type,model,link,ntrials,nsuccess,
               covariates.matrix.p,covariates.matrix.scalef,
               offset.p,offset.scalef) {
   nobs <- nrow(covariates.matrix.p) 
   output <- Model.Binary(parameter,model.type,model,link,ntrials,
                   covariates.matrix.p,covariates.matrix.scalef,
                   offset.p,offset.scalef) 
   probabilities <- output$probabilities 
# Calculation of log likelihood
      loglikelihood <- 0 
      for ( i in 1:nobs) { probability <- probabilities[[i]]
         probability <- probability*(probability>1.e-14) + 1.e-14*(probability<=1.e-14)
         vsuccess    <- nsuccess[[i]] 
         if (is.na(probability[1])==TRUE) { loglikelihood <- -1.e+20 
                 } else { loglikelihood  <- loglikelihood + 
                                t(log(probability))%*%vsuccess } } # end of for loop
   return(loglikelihood) }
