print.sightest <-
function(x,...){
    cat("Call:\n")
    print(x$call)
    cat("\n------------------- SIGHTABILITY MODEL ---------------------\n")
    print(x$sight.model)
    cat("\n----------------- Population Survey data  ----------------\n")
    cat("\n Stratum Sampling Information\n")
    print(x$samp)
   
    cat("\n Number of animals seen in each stratum\n")
    ua <- tapply(x$odat$total,x$odat$stratum, sum)
    print(format(ua, big.mark=","), quote = FALSE)
    Tot.seen <- sum(ua)
    cat("\n-------------- POPULATION ESTIMATE (",100*(1-x$alpha),"% CI) ----------------\n")
    z <- qnorm(1-x$alpha/2)
    if(x$CI.method == "normal"){
      temp <- rep(x$est[1], 3)+c(0, -z, z)*sqrt(x$est[2])
    }
    else{  #CI under lognormal assumption (see Wong p. 65-67)
  
  #Estimated number of animals that were not seen (assumed to be lognormally distributed)
    tau.m.T <- x$est[1]-Tot.seen 
    cv2 <- x$est[2]/(tau.m.T)^2
    cfact <- exp(z*sqrt(log(1+cv2)))
    temp <- rep(tau.m.T, 3)*c(1, (1/cfact)*sqrt(1+cv2), cfact*sqrt(1+cv2))+ rep(Tot.seen, 3)
    }
    names(temp) <- NULL
    temp <- format(round(temp, 0), big.mark=",")
    cat("\n")
    temp2 <- paste("tau.hat = ", temp[1], ";  ",100*(1-x$alpha), "% CI = (", temp[2], ", ", temp[3], ")")
    print(temp2[1], quote = FALSE)
    cat("\n")
    cat("\n------------------  SE(tau.hat) --------------------------------\n")
    cat("Variance method: "); print(x$var.method)
    SE <- sqrt(x$est[2])
    names(SE) <- "SE"
    print(format(round(SE, 0),  big.mark=","), quote = FALSE)
    cat("\n-------------- Variance Components -------------------\n")
    cat("\n")
    print(format(round(x$est[3:5], 0), big.mark=","), quote = FALSE)
}
