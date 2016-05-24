summary.sightest <-
function(object,...){
  # total number of animals seen
    Tot.seen <- sum(object$odat$total) 
  
  # z statistic
    z <- qnorm(1-object$alpha/2)
  
  # Normal assumption confidence interval  
    if(object$CI.method == "normal"){
      temp <- rep(object$est[1], 3)+c(0, -z, z)*sqrt(object$est[2])
    }
    else{  
 #CI under lognormal assumption (see Wong p. 65-67)
 #Estimated number of animals that were not seen (assumed to be lognormally distributed)
    tau.m.T <- object$est[1]-Tot.seen 
    cv2 <- object$est[2]/(tau.m.T)^2
    cfact <- exp(z*sqrt(log(1+cv2)))
    temp <- rep(tau.m.T, 3)*c(1,(1/cfact)*sqrt(1+cv2), cfact*sqrt(1+cv2))+ rep(Tot.seen, 3)
    }
    names(temp) <- NULL
    temp <- format(round(temp, 0), big.mark=",")
    cat("\n")
    temp2 <- paste("tau.hat = ", temp[1], ";  ",100*(1-object$alpha), "% CI = (", temp[2], ", ", temp[3], ")")
    print(format(temp2[1], big.mark=","), quote = FALSE)
    out.summary <- list(tau.hat = temp[1], lcl = temp[2], ucl = temp[3])  
}
