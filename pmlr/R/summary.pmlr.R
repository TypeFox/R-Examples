summary.pmlr <- function(object, ...) {
   print(object$call)
   if(object$method == "wald") cat("\nWald confidence intervals and p-values\n\n") else cat("\nLikelihood confidence intervals and p-values\n\n")
   p <- dim(object$coefficients)[2]
   J <- dim(object$coefficients)[3]
   output1 <- array(data = NA, dim = c(p,9,J))
   dimnames(output1) <- list(dimnames(object$coefficients)[[2]], c("coef","se(coef)","lower(coef)","upper(coef)","OR","lower(OR)","upper(OR)","Chi-Square","Pr > ChiSq"), dimnames(object$coefficients)[[3]])
   for (i in 1:J) {
      output1[,1,i] <- ifelse(is.infinite(object$separation[,,i]), object$separation[,,i], object$coefficients[,,i])
      output1[,2,i] <- ifelse(is.infinite(object$separation[,,i]), NA, sqrt(diag(object$var[,,i])))
      output1[,3,i] <- object$CI.lower[,,i]
      output1[,4,i] <- object$CI.upper[,,i]
      if (object$method =="wald") {
         output1[,3,i] <- ifelse(is.infinite(object$separation[,,i]), NA, object$CI.lower[,,i]) 
         output1[,4,i] <- ifelse(is.infinite(object$separation[,,i]), NA, object$CI.upper[,,i])
      }
      output1[,5,i] <- ifelse(is.infinite(object$separation[,,i]), exp(object$separation[,,i]), exp(object$coefficients[,,i]))
      output1[,6,i] <- exp(object$CI.lower[,,i])
      output1[,7,i] <- exp(object$CI.upper[,,i]) 
      if (object$method =="wald") {
         output1[,6,i] <- ifelse(is.infinite(object$separation[,,i]), NA, exp(object$CI.lower[,,i])) 
         output1[,7,i] <- ifelse(is.infinite(object$separation[,,i]), NA, exp(object$CI.upper[,,i]))
      }
      output1[,8,i] <- object$stat[,,i]        
      output1[,9,i] <- object$pval[,,i]
      if (object$method =="wald") {
         output1[,8,i] <- ifelse(is.infinite(object$separation[,,i]), NA, object$stat[,,i]) 
         output1[,9,i] <- ifelse(is.infinite(object$separation[,,i]), NA, object$pval[,,i])
      }
   }
   print(output1)
   if (object$joint) {
      output2 <- array(data = NA, dim = c(p,2,2))
      dimnames(output2) <- list(dimnames(object$coefficients)[[2]], c("Chi-Square","Pr > ChiSq"), c("H_0: b_{i,1} = b_{i,2} = ... = b_{i,J} = 0","H_0: b_{i,1} = b_{i,2} = ... = b_{i,J}"))      
      output2[,1,1] <- object$stat.joint[,,1]      
      output2[,2,1] <- object$pval.joint[,,1]      
      output2[,1,2] <- object$stat.joint[,,2]      
      output2[,2,2] <- object$pval.joint[,,2]      
      cat("H_0: b_{1i} = ... = b_{Ji} = 0\n\n")
      print(output2[,,1])
      cat("\nH_0: b_{1i} = ... = b_{Ji}\n\n")
      print(output2[,,2])      
   }   
}
                                                      