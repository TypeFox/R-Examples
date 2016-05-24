ICA.ContCont <- function(T0S0, T1S1, T0T0=1, T1T1=1, S0S0=1, S1S1=1, 
                         T0T1=seq(-1, 1, by=.1), T0S1=seq(-1, 1, by=.1), T1S0=seq(-1, 1, by=.1), 
                         S0S1=seq(-1, 1, by=.1)) { 
  
  T0S0_hier <- T0S0[1]
  T1S1_hier <- T1S1[1]
  
  Results <- na.exclude(matrix(NA, 1, 9))  
  colnames(Results) <- c("T0T1", "T0S0", "T0S1", "T1S0", "T1S1", "S0S1", "ICA", "Sigma.Delta.T", "delta") 
  combins <- expand.grid(T0T1, T0S0_hier, T0S1, T1S0, T1S1_hier, S0S1)
  lengte <- dim(combins)[1]
  
  if (length(T0S0)>1){
    if (length(T0S0)<lengte){stop("The specified vector for T0S0 should be larger than ", lengte) }
    T0S0_vector <- T0S0[1:lengte] # sample
    combins[,2] <- T0S0_vector
  }
  if (length(T1S1)>1){
    if (length(T1S1)<lengte){stop("The specified vector for T1S1 should be larger than ", lengte) }
        T1S1_vector <- T1S1[1:lengte] # sample
    combins[,5] <- T1S1_vector
  }
  
  
  for (i in 1: nrow(combins)) {   
    T0T1 <- combins[i, 1]  
    T0S0 <- combins[i, 2] 
    T0S1 <- combins[i, 3]  
    T1S0 <- combins[i, 4]  
    T1S1 <- combins[i, 5]  
    S0S1 <- combins[i, 6] 
    Sigma_c <- diag(4)         
    Sigma_c[2,1] <- Sigma_c[1,2] <- T0T1 * (sqrt(T0T0)*sqrt(T1T1))
    Sigma_c[3,1] <- Sigma_c[1,3] <- T0S0 * (sqrt(T0T0)*sqrt(S0S0))
    Sigma_c[4,1] <- Sigma_c[1,4] <- T0S1 * (sqrt(T0T0)*sqrt(S1S1))
    Sigma_c[3,2] <- Sigma_c[2,3] <- T1S0 * (sqrt(T1T1)*sqrt(S0S0))
    Sigma_c[4,2] <- Sigma_c[2,4] <- T1S1 * (sqrt(T1T1)*sqrt(S1S1))
    Sigma_c[4,3] <- Sigma_c[3,4] <- S0S1 * (sqrt(S0S0)*sqrt(S1S1))
    Sigma_c[1,1] <- T0T0
    Sigma_c[2,2] <- T1T1
    Sigma_c[3,3] <- S0S0
    Sigma_c[4,4] <- S1S1
    Cor_c <- cov2cor(Sigma_c)
    Min.Eigen.Cor <- try(min(eigen(Cor_c)$values), TRUE)  
    
    if (Min.Eigen.Cor > 0) {
      ICA <- ((sqrt(S0S0*T0T0)*Cor_c[3,1])+(sqrt(S1S1*T1T1)*Cor_c[4,2])-(sqrt(S0S0*T1T1)*Cor_c[3,2])-(sqrt(S1S1*T0T0)*Cor_c[4,1]))/(sqrt((T0T0+T1T1-(2*sqrt(T0T0*T1T1)*Cor_c[2,1]))*(S0S0+S1S1-(2*sqrt(S0S0*S1S1)*Cor_c[4,3]))))
      if ((is.finite(ICA))==TRUE){
        sigma.delta.T <- T0T0 + T1T1 - (2 * sqrt(T0T0*T1T1) * Cor_c[2,1])
        delta <- sigma.delta.T * (1-(ICA**2))
        results.part <- as.vector(cbind(T0T1, T0S0, T0S1, T1S0, T1S1, S0S1, ICA, sigma.delta.T, delta))
        Results <- rbind(Results, results.part)
        rownames(Results) <- NULL}
    }
  }
  Results <- data.frame(Results)
  rownames(Results) <- NULL
  Total.Num.Matrices <- nrow(combins)
  
  fit <- 
    list(Total.Num.Matrices=Total.Num.Matrices, Pos.Def=Results[,1:6], ICA=Results$ICA, GoodSurr=Results[,7:9], Call=match.call())
  
  class(fit) <- "ICA.ContCont"
  fit
  
}


plot.ICA.ContCont <- function(x, Xlab.ICA, Main.ICA, Type="Percent", Labels=FALSE, ICA=TRUE, Good.Surr=FALSE, Main.Good.Surr, 
                              Par=par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1)), col, ...){   

 
    Object <- x 
    if (missing(Xlab.ICA)) {Xlab.ICA <- expression(rho[Delta])}
    if (missing(Main.ICA)) {Main.ICA="ICA"} 
    if (missing(col)) {col <- c(8)}
    
    if (ICA==TRUE){
    
    dev.new() 
    par=Par  
    if (Type=="Freq"){
        h <- hist(Object$ICA, ...)
        h$density <- h$counts/sum(h$counts)
        cumulMidPoint <- ecdf(x=Object$ICA)(h$mids)
        labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
        
        if (Labels==FALSE){
          plot(h,freq=T, xlab=Xlab.ICA, ylab="Frequency", col=col, main=Main.ICA)
           }
        if (Labels==TRUE){
          plot(h,freq=T, xlab=Xlab.ICA, ylab="Frequency", col=col, main=Main.ICA, labels=labs)
           }
        }
    
    if (Type=="Percent"){
        h <- hist(Object$ICA, ...)
        h$density <- h$counts/sum(h$counts)
        cumulMidPoint <- ecdf(x=Object$ICA)(h$mids)
        labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
        
        if (Labels==FALSE){
          plot(h,freq=F, xlab=Xlab.ICA, ylab="Percentage", col=col, main=Main.ICA)
            }
        if (Labels==TRUE){
          plot(h,freq=F, xlab=Xlab.ICA, ylab="Percentage", col=col, main=Main.ICA, labels=labs)
            }
        }
    
    if (Type=="CumPerc"){
      h <- hist(Object$ICA, breaks=length(Object$ICA), ...)
      h$density <- h$counts/sum(h$counts)
      cumulative <- cumsum(h$density)
      plot(x=h$mids, y=cumulative, xlab=Xlab.ICA, ylab="Cumulative percentage", col=0, main=Main.ICA)
      lines(x=h$mids, y=cumulative)
         }
    }
    
    
    if (Good.Surr==TRUE){
    
    if (missing(Main.Good.Surr)) {Main.Good.Surr = " "}  
    par=Par
    if (Type=="Freq"){
      h <- hist(Object$GoodSurr$delta, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$GoodSurr$delta)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h,freq=T, xlab=expression(delta), ylab="Frequency", main=Main.Good.Surr, col=col)
      }
      if (Labels==TRUE){
        plot(h,freq=T, xlab=expression(delta), ylab="Frequency", col=col, labels=labs, main=Main.Good.Surr)
      }
    }
    
    if (Type=="Percent"){
      h <- hist(Object$GoodSurr$delta, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$GoodSurr$delta)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h,freq=F, xlab=expression(delta), ylab="Percentage", col=col, main=Main.Good.Surr)
      }
      if (Labels==TRUE){
        plot(h,freq=F, xlab=expression(delta), ylab="Percentage", col=col, labels=labs, main=Main.Good.Surr)
      }
    }
    
    if (Type=="CumPerc"){
      h <- hist(Object$GoodSurr$delta, breaks=length(Object$GoodSurr$delta), ...)
      h$density <- h$counts/sum(h$counts)
      cumulative <- cumsum(h$density)
      plot(x=h$mids, y=cumulative, xlab=expression(delta), ylab="Cumulative percentage", col=0, main=Main.Good.Surr)
      lines(x=h$mids, y=cumulative)
    }
    
  
    }    

}

summary.ICA.ContCont <- function(object, ..., Object){
  
  
  if (missing(Object)){Object <- object} 
  
  mode <- function(data) {
    x <- data
    z <- density(x)
    mode_val <- z$x[which.max(z$y)]
    fit <- list(mode_val= mode_val)
  }
  
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n\n# Total number of matrices that can be formed by the specified vectors and/or scalars")
  cat("\n# of the correlations in the function call")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  cat(Object$Total.Num.Matrices)
  cat("\n\n# Total number of positive definite matrices")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  cat(nrow(Object$Pos.Def))
  cat("\n\n\n# Causal-inference (ICA) results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")  
  cat("Mean (SD) ICA: ", format(round(mean(Object$ICA), 4), nsmall = 4), " (", format(round(sd(Object$ICA), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(Object$ICA), 4), nsmall = 4), "; max: ",  format(round(max(Object$ICA), 4), nsmall = 4), "]", sep="")
  cat("\nMode ICA: ", format(round(mode(Object$ICA)$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the ICA distribution: \n\n")
  quant <- quantile(Object$ICA, probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)
}

