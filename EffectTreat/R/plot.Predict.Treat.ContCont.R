plot.Predict.Treat.ContCont <- plot.Predict.Treat.Multivar.ContCont <- function(x, Xlab, Main, Mean.T0T1=FALSE, Median.T0T1=TRUE,  
    Specific.T0T1="none", alpha=0.05, Cex.Legend=1, ...){ 

  Object <- x 

  if (class(Object)=="Predict.Treat.ContCont"){
  
  if (missing(Xlab)) {Xlab <- expression(paste(Delta, "T"[j], "|S"[j]))}
  if (missing(Main)) {Main=" "} 
  
  min_T0T1 <- min(Object$T0T1)
  max_T0T1 <- max(Object$T0T1)
  
  if (Specific.T0T1!="none"){
    if ((Specific.T0T1 < min_T0T1) | (Specific.T0T1 > max_T0T1)){
      stop("The requested Specific.T0T1 value lies outside the range for which valid solutions were obtained. Please specify a value >=", min_T0T1, " and <= ", max_T0T1)
        }
    }
    
  Var_Delta.T_S <- Object$Var_Delta.T_S
  mean_Pred_T <- Object$Pred_T[1]
  crit_val <- qnorm(c(alpha/2), mean=0, sd=1, lower.tail=FALSE)
  
  if (class(Object)=="Predict.Treat.ContCont") {
   
      min_rho_T0T1 <- min(Object$T0T1)
      max_rho_T0T1 <- max(Object$T0T1)
      median_rho_T0T1 <- median(Object$T0T1)
      mean_rho_T0T1 <- mean(Object$T0T1)
      
      min_Var_Delta.T <- Object$Var_Delta.T[Object$T0T1==min_rho_T0T1]
      max_Var_Delta.T <- Object$Var_Delta.T[Object$T0T1==max_rho_T0T1]
      median_Var_Delta.T <- Object$Var_Delta.T[which.min(abs(Object$T0T1-median_rho_T0T1))]
      mean_Var_Delta.T <- Object$Var_Delta.T[which.min(abs(Object$T0T1-mean_rho_T0T1))]
      
      min_SD_Delta.T_givenS <- sqrt(min_Var_Delta.T*(1-((Object$PCA[Object$T0T1==min_rho_T0T1])**2)))
      max_SD_Delta.T_givenS <- sqrt(max_Var_Delta.T*(1-((Object$PCA[Object$T0T1==max_rho_T0T1])**2)))
      median_SD_Delta.T_givenS <- sqrt(median_Var_Delta.T*(1-((Object$PCA[which.min(abs(Object$T0T1-median_rho_T0T1))])**2)))
      mean_SD_Delta.T_givenS <- sqrt(mean_Var_Delta.T*(1-((Object$PCA[which.min(abs(Object$T0T1-mean_rho_T0T1))])**2)))
      rm(min_Var_Delta.T, max_Var_Delta.T, median_Var_Delta.T, mean_Var_Delta.T)

      x <- ((seq(-4,4,length=1000)*mean_SD_Delta.T_givenS))*1 + mean_Pred_T
      hx <- dnorm(x, mean_Pred_T, mean_SD_Delta.T_givenS)
      plot(x, hx*crit_val, type="l", xlab=Xlab, ylab="", main=Main, lwd=2, col=0, ...)  
            
      # mean rho_T0T1
      if (Mean.T0T1==TRUE){
        x <- ((seq(-4,4,length=1000)*mean_SD_Delta.T_givenS))*1 + mean_Pred_T
        hx <- dnorm(x, mean_Pred_T, mean_SD_Delta.T_givenS)
        lines(x=x, y=hx, col=2, lwd=2)
        
        SD_hier <- mean_SD_Delta.T_givenS
        col_hier <- 2
        x1 <- mean_Pred_T - (SD_hier * crit_val)
        x2 <- mean_Pred_T + (SD_hier * crit_val)
        sam <- cbind(x, hx)
        y1 <- sam[,2][which.min(abs(sam[,1]-x1))]
        y2 <- sam[,2][which.min(abs(sam[,1]-x2))]
        
        segments(x0 = x1, y0 = 0, x1 = x1, y1 = y1, col=col_hier, lwd=2, lty=2)
        segments(x0 = x2, y0 = 0, x1 = x2, y1 = y2, col=col_hier, lwd=2, lty=2)
      }    
      
      # median rho_T0T1
      if (Median.T0T1==TRUE){
        x <- ((seq(-4,4,length=1000)*median_SD_Delta.T_givenS))*1 + mean_Pred_T
        hx <- dnorm(x, mean_Pred_T, median_SD_Delta.T_givenS)
        lines(x=x, y=hx, col=3, lwd=2)

        SD_hier <- median_SD_Delta.T_givenS
        col_hier <- 3
        x1 <- mean_Pred_T - (SD_hier * crit_val)
        x2 <- mean_Pred_T + (SD_hier * crit_val)
        sam <- cbind(x, hx)
        y1 <- sam[,2][which.min(abs(sam[,1]-x1))]
        y2 <- sam[,2][which.min(abs(sam[,1]-x2))]
        
        segments(x0 = x1, y0 = 0, x1 = x1, y1 = y1, col=col_hier, lwd=2, lty=2)
        segments(x0 = x2, y0 = 0, x1 = x2, y1 = y2, col=col_hier, lwd=2, lty=2)
        
      }
      
      # specific value for rho_T0T1
      if (Specific.T0T1!="none"){
       
        user_req_Var_Delta.T_givenS <- Object$Var_Delta.T_S[which.min(abs(Object$T0T1-Specific.T0T1))]
        PCA_val <- Object$PCA[which.min(abs(Object$T0T1-Specific.T0T1))]
        
        user_req_SD_Delta.T_givenS <- sqrt(user_req_Var_Delta.T_givenS)
        
        x <- ((seq(-4,4,length=1000)*user_req_SD_Delta.T_givenS))*1 + mean_Pred_T
        hx <- dnorm(x, mean_Pred_T, user_req_SD_Delta.T_givenS)
        lines(x=x, y=hx, col=4, lwd=2)

        SD_hier <- user_req_SD_Delta.T_givenS
        col_hier <- 4
        x1 <- mean_Pred_T - (SD_hier * crit_val)
        x2 <- mean_Pred_T + (SD_hier * crit_val)
        sam <- cbind(x, hx)
        y1 <- sam[,2][which.min(abs(sam[,1]-x1))]
        y2 <- sam[,2][which.min(abs(sam[,1]-x2))]
        
        segments(x0 = x1, y0 = 0, x1 = x1, y1 = y1, col=col_hier, lwd=2, lty=2)
        segments(x0 = x2, y0 = 0, x1 = x2, y1 = y2, col=col_hier, lwd=2, lty=2)
        
        
        
      }
      
      abline(v = mean_Pred_T, lty=2)
      
      
      if (Mean.T0T1==TRUE & Median.T0T1==TRUE & Specific.T0T1!="none"){
      legend("topright", inset=.05, legend=c(bquote(paste(rho[T0T1], "="~.(mean_rho_T0T1))), 
      bquote(paste(rho[T0T1], "="~.(median_rho_T0T1))), bquote(paste(rho[T0T1], "="~.(Specific.T0T1))), expression()), 
             lwd=2, col=c(2, 3, 4), cex=Cex.Legend)
      }
      
      if (Mean.T0T1==FALSE & Median.T0T1==TRUE & Specific.T0T1!="none"){
      legend("topright", inset=.05, legend=c( 
      bquote(paste(rho[T0T1], "="~.(median_rho_T0T1))), bquote(paste(rho[T0T1], "="~.(Specific.T0T1))), expression()), 
               lwd=2, col=c(3, 4), cex=Cex.Legend)
      }
      
      if (Mean.T0T1==TRUE & Median.T0T1==TRUE & Specific.T0T1=="none"){
        legend("topright", inset=.05, legend=c(bquote(paste(rho[T0T1], "="~.(mean_rho_T0T1))), 
                                               bquote(paste(rho[T0T1], "="~.(median_rho_T0T1))), expression()), 
               lwd=2, col=c(2, 3), cex=Cex.Legend)
      }

      if (Mean.T0T1==FALSE & Median.T0T1==TRUE & Specific.T0T1=="none"){
        legend("topright", inset=.05, legend=c( 
        bquote(paste(rho[T0T1], "="~.(median_rho_T0T1))), expression()), 
               lwd=2, col=c(3), cex=Cex.Legend)
      }
      
      if (Mean.T0T1==FALSE & Median.T0T1==FALSE & Specific.T0T1!="none"){
        legend("topright", inset=.05, legend=c(bquote(paste(rho[T0T1], "="~.(Specific.T0T1))), expression()), 
               lwd=2, col=c(4), cex=Cex.Legend)
      }
      
      if (Mean.T0T1==TRUE & Median.T0T1==FALSE & Specific.T0T1=="none"){
        legend("topright", inset=.05, legend=c(bquote(paste(rho[T0T1], "="~.(mean_rho_T0T1))), expression()), 
               lwd=2, col=c(2), cex=Cex.Legend)
      }

      if (Mean.T0T1==TRUE & Median.T0T1==FALSE & Specific.T0T1!="none"){
        legend("topright", inset=.05, legend=c(bquote(paste(rho[T0T1], "="~.(mean_rho_T0T1))), 
                                               bquote(paste(rho[T0T1], "="~.(Specific.T0T1))), expression()), 
               lwd=2, col=c(2, 4), cex=Cex.Legend)
      }
      
      
   }
  }
  
  
  if (class(Object)=="Predict.Treat.Multivar.ContCont"){
    
    if (missing(Xlab)) {Xlab <- expression(paste(Delta, "T"[j], "|S"[j]))}
    if (missing(Main)) {Main=" "} 
    
    min_T0T1 <- min(Object$T0T1)
    max_T0T1 <- max(Object$T0T1)
    
    if (Specific.T0T1!="none"){
      if ((Specific.T0T1 < min_T0T1) | (Specific.T0T1 > max_T0T1)){
        stop("The requested Specific.T0T1 value lies outside the range for which valid solutions were obtained. Please specify a value >=", min_T0T1, " and <= ", max_T0T1)
      }
    }
    
    Var_Delta.T_S <- Object$Var_Delta.T_S
    mean_Pred_T <- Object$Pred_T[1]
    crit_val <- qnorm(c(alpha/2), mean=0, sd=1, lower.tail=FALSE)
    
    if (class(Object)=="Predict.Treat.Multivar.ContCont") {
      
      min_rho_T0T1 <- min(Object$T0T1)
      max_rho_T0T1 <- max(Object$T0T1)
      median_rho_T0T1 <- median(Object$T0T1)
      mean_rho_T0T1 <- mean(Object$T0T1)
      
      min_SD_Delta.T_givenS <- sqrt(Object$Var_Delta.T_S[Object$T0T1==min_rho_T0T1])
      max_SD_Delta.T_givenS <- sqrt(Object$Var_Delta.T_S[Object$T0T1==max_rho_T0T1])
      median_SD_Delta.T_givenS <- sqrt(Object$Var_Delta.T_S[which.min(abs(Object$T0T1-median_rho_T0T1))])
      mean_SD_Delta.T_givenS <- sqrt(Object$Var_Delta.T_S[which.min(abs(Object$T0T1-mean_rho_T0T1))])
      
      x <- ((seq(-4,4,length=1000)*mean_SD_Delta.T_givenS))*1 + mean_Pred_T
      hx <- dnorm(x, mean_Pred_T, mean_SD_Delta.T_givenS)
      plot(x, hx*crit_val, type="l", xlab=Xlab, ylab="", main=Main, lwd=2, col=0, ...)  
      
      # mean rho_T0T1
      if (Mean.T0T1==TRUE){
        x <- ((seq(-4,4,length=1000)*mean_SD_Delta.T_givenS))*1 + mean_Pred_T
        hx <- dnorm(x, mean_Pred_T, mean_SD_Delta.T_givenS)
        lines(x=x, y=hx, col=2, lwd=2)
        
        SD_hier <- mean_SD_Delta.T_givenS
        col_hier <- 2
        x1 <- mean_Pred_T - (SD_hier * crit_val)
        x2 <- mean_Pred_T + (SD_hier * crit_val)
        sam <- cbind(x, hx)
        y1 <- sam[,2][which.min(abs(sam[,1]-x1))]
        y2 <- sam[,2][which.min(abs(sam[,1]-x2))]
        
        segments(x0 = x1, y0 = 0, x1 = x1, y1 = y1, col=col_hier, lwd=2, lty=2)
        segments(x0 = x2, y0 = 0, x1 = x2, y1 = y2, col=col_hier, lwd=2, lty=2)
      }    
      
      # median rho_T0T1
      if (Median.T0T1==TRUE){
        x <- ((seq(-4,4,length=1000)*median_SD_Delta.T_givenS))*1 + mean_Pred_T
        hx <- dnorm(x, mean_Pred_T, median_SD_Delta.T_givenS)
        lines(x=x, y=hx, col=3, lwd=2)
        
        SD_hier <- median_SD_Delta.T_givenS
        col_hier <- 3
        x1 <- mean_Pred_T - (SD_hier * crit_val)
        x2 <- mean_Pred_T + (SD_hier * crit_val)
        sam <- cbind(x, hx)
        y1 <- sam[,2][which.min(abs(sam[,1]-x1))]
        y2 <- sam[,2][which.min(abs(sam[,1]-x2))]
        
        segments(x0 = x1, y0 = 0, x1 = x1, y1 = y1, col=col_hier, lwd=2, lty=2)
        segments(x0 = x2, y0 = 0, x1 = x2, y1 = y2, col=col_hier, lwd=2, lty=2)
        
      }
      
      # specific value for rho_T0T1
      if (Specific.T0T1!="none"){
        
        user_req_Var_Delta.T_givenS <- Object$Var_Delta.T_S[which.min(abs(Object$T0T1-Specific.T0T1))]
        PCA_val <- Object$PCA[which.min(abs(Object$T0T1-Specific.T0T1))]
        
        user_req_SD_Delta.T_givenS <- sqrt(user_req_Var_Delta.T_givenS)
        
        x <- ((seq(-4,4,length=1000)*user_req_SD_Delta.T_givenS))*1 + mean_Pred_T
        hx <- dnorm(x, mean_Pred_T, user_req_SD_Delta.T_givenS)
        lines(x=x, y=hx, col=4, lwd=2)
        
        SD_hier <- user_req_SD_Delta.T_givenS
        col_hier <- 4
        x1 <- mean_Pred_T - (SD_hier * crit_val)
        x2 <- mean_Pred_T + (SD_hier * crit_val)
        sam <- cbind(x, hx)
        y1 <- sam[,2][which.min(abs(sam[,1]-x1))]
        y2 <- sam[,2][which.min(abs(sam[,1]-x2))]
        
        segments(x0 = x1, y0 = 0, x1 = x1, y1 = y1, col=col_hier, lwd=2, lty=2)
        segments(x0 = x2, y0 = 0, x1 = x2, y1 = y2, col=col_hier, lwd=2, lty=2)
        
      }
      
      abline(v = mean_Pred_T, lty=2)
      
      
      if (Mean.T0T1==TRUE & Median.T0T1==TRUE & Specific.T0T1!="none"){
        legend("topright", inset=.05, legend=c(bquote(paste(rho[T0T1], "="~.(mean_rho_T0T1))), 
                                               bquote(paste(rho[T0T1], "="~.(median_rho_T0T1))), bquote(paste(rho[T0T1], "="~.(Specific.T0T1))), expression()), 
               lwd=2, col=c(2, 3, 4), cex=Cex.Legend)
      }
      
      if (Mean.T0T1==FALSE & Median.T0T1==TRUE & Specific.T0T1!="none"){
        legend("topright", inset=.05, legend=c( 
          bquote(paste(rho[T0T1], "="~.(median_rho_T0T1))), bquote(paste(rho[T0T1], "="~.(Specific.T0T1))), expression()), 
          lwd=2, col=c(3, 4), cex=Cex.Legend)
      }
      
      if (Mean.T0T1==TRUE & Median.T0T1==TRUE & Specific.T0T1=="none"){
        legend("topright", inset=.05, legend=c(bquote(paste(rho[T0T1], "="~.(mean_rho_T0T1))), 
                                               bquote(paste(rho[T0T1], "="~.(median_rho_T0T1))), expression()), 
               lwd=2, col=c(2, 3), cex=Cex.Legend)
      }
      
      if (Mean.T0T1==FALSE & Median.T0T1==TRUE & Specific.T0T1=="none"){
        legend("topright", inset=.05, legend=c( 
          bquote(paste(rho[T0T1], "="~.(median_rho_T0T1))), expression()), 
          lwd=2, col=c(3), cex=Cex.Legend)
      }
      
      if (Mean.T0T1==FALSE & Median.T0T1==FALSE & Specific.T0T1!="none"){
        legend("topright", inset=.05, legend=c(bquote(paste(rho[T0T1], "="~.(Specific.T0T1))), expression()), 
               lwd=2, col=c(4), cex=Cex.Legend)
      }
      
      if (Mean.T0T1==TRUE & Median.T0T1==FALSE & Specific.T0T1=="none"){
        legend("topright", inset=.05, legend=c(bquote(paste(rho[T0T1], "="~.(mean_rho_T0T1))), expression()), 
               lwd=2, col=c(2), cex=Cex.Legend)
      }
      
      if (Mean.T0T1==TRUE & Median.T0T1==FALSE & Specific.T0T1!="none"){
        legend("topright", inset=.05, legend=c(bquote(paste(rho[T0T1], "="~.(mean_rho_T0T1))), 
                                               bquote(paste(rho[T0T1], "="~.(Specific.T0T1))), expression()), 
               lwd=2, col=c(2, 4), cex=Cex.Legend)
      }
      
      
    }
  }
  
  
  
}