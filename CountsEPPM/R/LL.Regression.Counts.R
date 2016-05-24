LL.Regression.Counts <-
function(parameter,model.type,model,list.counts,
               covariates.matrix.mean,covariates.matrix.variance,
               offset.mean,offset.variance,ltvalue,utvalue,
               scale.factor.model,fixed.b) {

   nobs   <- nrow(covariates.matrix.mean) 
   vnmax  <- sapply(list.counts,length) - 1
   output <- Model.Counts(parameter,model.type,model,covariates.matrix.mean,
                   covariates.matrix.variance,offset.mean,offset.variance,
                   scale.factor.model,fixed.b,vnmax) 

# 1.e-323 = exp(-743.7469) is the smallest probability value not to give -Inf
# 1.e-325 has been set as value for -infinity for log likelihoods
      small.loglikelihood    <- log(1.e-323)
      infinity.loglikelihood <- log(1.e-325)

# Calculation of log likelihood
      loglikelihood <- 0 
      for ( i in 1:nobs) { probability <- output$probabilities[[i]]
         count <- list.counts[[i]] 
         nmax1 <- vnmax[i] + 1

# With overly large values of the parameters and small counts it is possible to have
# probabilities exactly 0 returned. Also LRTruncation sets the probabilities below
# ltvalue and above utvalue to NA.

# rounding error can cause the sum of the probabilities to <0 or >1
# so rounded to 10 decimal places
        wks <- round(sum(probability),digits=10)
        if ((is.finite(wks)==FALSE) | (wks<=0) | (wks>1)) {
            wk.loglikelihood <- sum(count)*small.loglikelihood  

                                                    } else {

           if ((is.na(ltvalue)==FALSE) | (is.na(utvalue)==FALSE)) { 
             if (is.na(ltvalue)==FALSE) { wks1 <- ltvalue + 2 
                                        } else { wks1 <- 1 }
             if (is.na(utvalue)==FALSE) { wks2 <- utvalue  
                                        } else { wks2 <- nmax1 }
             rev.probability <- LRTruncation(probability,ltvalue,utvalue) 
             wks  <- round(sum(rev.probability[wks1:wks2]),digits=10)
             if (is.finite(wks)==FALSE) { wks <- 0 }
             wkv2 <- count[wks1:wks2]
             if ((is.na(utvalue)==FALSE) & (wks!=1)) { 
                wk.loglikelihood  <- infinity.loglikelihood 
                                                     } else {
                wkv1 <- log(rev.probability[wks1:wks2])
                nlen <- length(wkv1)
                wkv1 <- sapply(1:nlen, function(j) 
                   if (is.finite(wkv1[j])==FALSE) { wkv1[j] <- small.loglikelihood 
                                                  } else { wkv1[j] <- wkv1[j] } )
                       wk.loglikelihood <- t(wkv1)%*%wkv2 } # end wks!=1
                                                                    } else {
                wkv1 <- log(probability)
                wkv1 <- sapply(1:nmax1, function(j) 
                   if (is.finite(wkv1[j])==FALSE) { wkv1[j] <- small.loglikelihood 
                                                  } else { wkv1[j] <- wkv1[j] } )
                wk.loglikelihood <- t(wkv1)%*%count 
                                            } # end of is.na(ltvalue) | is.na(utvalue)
                                            } # end of if ((is.finite(wks)==FALSE) | (wks<=0) | (wks>1)) 
                loglikelihood  <- loglikelihood + wk.loglikelihood 
                                            } # end of for loop

   if (is.finite(loglikelihood)==FALSE) { 
        loglikelihood <-  sum(vnmax)*infinity.loglikelihood }

   return(loglikelihood)                          }
