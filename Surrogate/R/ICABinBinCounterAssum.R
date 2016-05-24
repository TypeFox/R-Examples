ICA.BinBin.CounterAssum <- function(x, r2_h_S0S1_min, r2_h_S0S1_max, r2_h_S0T1_min, r2_h_S0T1_max, 
                                    r2_h_T0T1_min, r2_h_T0T1_max, r2_h_T0S1_min, r2_h_T0S1_max, Monotonicity="General",
                                    Type="Freq", MainPlot=" ", Cex.Legend=1, Cex.Position="topright", ...) {
  
  if (class(x)!="ICA.BinBin") {stop("The function R2HBinBinCounterAssumpt should be applied to an object of 
                                    class ICA.BinBin.")}
  
  sub <- data.frame(cbind(x$Pi.Vectors, x$R2_H, x$Theta_T, x$Theta_S))
  
  if ((Monotonicity=="No" | Monotonicity=="True.Endp" | Monotonicity=="Surr.Endp" | Monotonicity=="Surr.True.Endp" | Monotonicity=="General")==FALSE){
    stop("The Monotonicity=... argument is not correctly specified \n")
  }
  
  if (Monotonicity=="No"){sub <- sub[sub$Monotonicity=="No",]}
  if (Monotonicity=="True.Endp"){sub <- sub[sub$Monotonicity=="True",]}
  if (Monotonicity=="Surr.Endp"){sub <- sub[sub$Monotonicity=="Surr",]}
  if (Monotonicity=="Surr.True.Endp"){sub <- sub[sub$Monotonicity=="SurrTrue",]}
  
  if (dim(sub)[1]==0) {stop("There are no valid observations for the specified ranges [min, max].")}
  
  #S0 S1
  cell_00 <- sub$Pi_0000 + sub$Pi_0100 + sub$Pi_1000 + sub$Pi_1100     
  cell_10 <- sub$Pi_0010 + sub$Pi_1010 + sub$Pi_1110 + sub$Pi_0110      
  cell_01 <- sub$Pi_0001 + sub$Pi_0101 + sub$Pi_1001 + sub$Pi_1101      
  cell_11 <- sub$Pi_1011 + sub$Pi_1111 + sub$Pi_0011 + sub$Pi_0111 
  sum_0 <- cell_00 + cell_10  
  sum_1 <- cell_01 + cell_11
  sum0_ <- cell_00 + cell_01 
  sum1_ <- cell_10 + cell_11
  
  R2_H_S0S1 <- matrix(NA, ncol = length(sub$Monotonicity)) 
  for (i in 1: length(sub$Monotonicity)){
    if (sub$Monotonicity[i]=="No" | sub$Monotonicity[i]=="True"){
      I <- (cell_00[i] * log2(cell_00[i]/(sum0_[i] * sum_0[i])))+
        (cell_10[i] * log2(cell_10[i]/(sum1_[i] * sum_0[i])))+  #
        (cell_01[i] * log2(cell_01[i]/(sum0_[i] * sum_1[i])))+
        (cell_11[i] * log2(cell_11[i]/(sum1_[i] * sum_1[i])))
      H_S <- - ((sum_0[i] * log2(sum_0[i])) + (sum_1[i] * log2(sum_1[i])))
      H_T <- - ((sum0_[i] * log2(sum0_[i])) + (sum1_[i] * log2(sum1_[i]))) 
      R2_H_S0S1[i] <- 
        as.numeric(I/(min(H_S, H_T))) 
    }
    
    if (sub$Monotonicity[i]=="Surr" | sub$Monotonicity[i]=="SurrTrue"){
      I <- (cell_00[i] * log2(cell_00[i]/(sum0_[i] * sum_0[i])))+
        (cell_01[i] * log2(cell_01[i]/(sum0_[i] * sum_1[i])))+
        (cell_11[i] * log2(cell_11[i]/(sum1_[i] * sum_1[i])))
      H_S <- - ((sum_0[i] * log2(sum_0[i])) + (sum_1[i] * log2(sum_1[i])))
      H_T <- - ((sum0_[i] * log2(sum0_[i])) + (sum1_[i] * log2(sum1_[i]))) 
      R2_H_S0S1[i] <- 
        as.numeric(I/(min(H_S, H_T))) 
    }  
  }
  
  #T0 T1
  cell_00 <- sub$Pi_0000 + sub$Pi_0010 + sub$Pi_0001 + sub$Pi_0011          
  cell_10 <- sub$Pi_1000 + sub$Pi_1010 + sub$Pi_1001 + sub$Pi_1011     
  cell_01 <- sub$Pi_0100 + sub$Pi_0101 + sub$Pi_0110 + sub$Pi_0111         
  cell_11 <- sub$Pi_1110 + sub$Pi_1101 + sub$Pi_1111 + sub$Pi_1100  
  sum_0 <- cell_00 + cell_10  
  sum_1 <- cell_01 + cell_11
  sum0_ <- cell_00 + cell_01 
  sum1_ <- cell_10 + cell_11
  
  R2_H_T0T1 <- matrix(NA, ncol = length(sub$Monotonicity))  
  
  for (i in 1: length(sub$Monotonicity)){
    
    if (sub$Monotonicity[i]=="No" | sub$Monotonicity[i]=="Surr"){
      
      I <- (cell_00[i] * log2(cell_00[i]/(sum0_[i] * sum_0[i])))+
        (cell_10[i] * log2(cell_10[i]/(sum1_[i] * sum_0[i])))+
        (cell_01[i] * log2(cell_01[i]/(sum0_[i] * sum_1[i])))+
        (cell_11[i] * log2(cell_11[i]/(sum1_[i] * sum_1[i])))  
      H_S <- - ((sum_0[i] * log2(sum_0[i])) + (sum_1[i] * log2(sum_1[i])))
      H_T <- - ((sum0_[i] * log2(sum0_[i])) + (sum1_[i] * log2(sum1_[i]))) 
      R2_H_T0T1[i] <- 
        I/(min(H_S, H_T)) 
    }
    
    if (sub$Monotonicity[i]=="True" | sub$Monotonicity[i]=="SurrTrue"){
      
      I <- (cell_00[i] * log2(cell_00[i]/(sum0_[i] * sum_0[i])))+
        (cell_01[i] * log2(cell_01[i]/(sum0_[i] * sum_1[i])))+
        (cell_11[i] * log2(cell_11[i]/(sum1_[i] * sum_1[i])))  
      H_S <- - ((sum_0[i] * log2(sum_0[i])) + (sum_1[i] * log2(sum_1[i])))
      H_T <- - ((sum0_[i] * log2(sum0_[i])) + (sum1_[i] * log2(sum1_[i]))) 
      R2_H_T0T1[i] <- 
        I/(min(H_S, H_T)) 
    }
    
  }
  
  
  #S0 T1
  cell_00 <- sub$Pi_0000 + sub$Pi_0001 + sub$Pi_1000 + sub$Pi_1001         
  cell_10 <- sub$Pi_0010 + sub$Pi_1010 + sub$Pi_1011 + sub$Pi_0011     
  cell_01 <- sub$Pi_0100 + sub$Pi_0101 + sub$Pi_1101 + sub$Pi_1100          
  cell_11 <- sub$Pi_1110 + sub$Pi_1111 + sub$Pi_0110 + sub$Pi_0111   
  sum_0 <- cell_00 + cell_10  
  sum_1 <- cell_01 + cell_11
  sum0_ <- cell_00 + cell_01 
  sum1_ <- cell_10 + cell_11
  
  I <- (cell_00 * log2(cell_00/(sum0_ * sum_0)))+
    (cell_10 * log2(cell_10/(sum1_ * sum_0)))+
    (cell_01 * log2(cell_01/(sum0_ * sum_1)))+
    (cell_11 * log2(cell_11/(sum1_ * sum_1)))
  H_S <- - ((sum_0 * log2(sum_0)) + (sum_1 * log2(sum_1)))
  H_T <- - ((sum0_ * log2(sum0_)) + (sum1_ * log2(sum1_))) 
  R2_H_S0T1 <- 
    I/(min(H_S, H_T)) 
  
  
  #S1 T0
  cell_00 <- sub$Pi_0000 + sub$Pi_0100 + sub$Pi_0010 + sub$Pi_0110          
  cell_10 <- sub$Pi_0001 + sub$Pi_0101 + sub$Pi_0011 + sub$Pi_0111      
  cell_01 <- sub$Pi_1000 + sub$Pi_1010 + sub$Pi_1110 + sub$Pi_1100         
  cell_11 <- sub$Pi_1001 + sub$Pi_1101 + sub$Pi_1011 + sub$Pi_1111
  sum_0 <- cell_00 + cell_10  
  sum_1 <- cell_01 + cell_11
  sum0_ <- cell_00 + cell_01 
  sum1_ <- cell_10 + cell_11
  
  I <- (cell_00 * log2(cell_00/(sum0_ * sum_0)))+
    (cell_10 * log2(cell_10/(sum1_ * sum_0)))+
    (cell_01 * log2(cell_01/(sum0_ * sum_1)))+
    (cell_11 * log2(cell_11/(sum1_ * sum_1)))
  H_S <- - ((sum_0 * log2(sum_0)) + (sum_1 * log2(sum_1)))
  H_T <- - ((sum0_ * log2(sum0_)) + (sum1_ * log2(sum1_))) 
  R2_H_S1T0 <- 
    I/(min(H_S, H_T)) 
  
  
  compA <- sub$x.R2_H
  compB <- as.numeric(R2_H_S0S1)
  compC <- as.numeric(R2_H_T0T1)
  compD <- as.numeric(R2_H_S0T1)
  compE <- as.numeric(R2_H_S1T0)
  compF <- sub$Monotonicity
  
  results <-
     data.frame(compA, compB, compC, compD, compE, as.character(compF))
  
  colnames(results) <- c("R2_H", "r2_h_S0S1", "r2_h_T0T1", "r2_h_S0T1", "r2_h_S1T0", "Monotonicity")
  
  subset <- results <- data.frame(results)
  
  subset <- subset[subset$r2_h_S0S1 > r2_h_S0S1_min & subset$r2_h_S0S1 < r2_h_S0S1_max,]
  subset <- subset[subset$r2_h_T0T1 > r2_h_T0T1_min & subset$r2_h_T0T1 < r2_h_T0T1_max,]
  subset <- subset[subset$r2_h_S0T1 > r2_h_S0T1_min & subset$r2_h_S0T1 < r2_h_S0T1_max,]
  subset <- subset[subset$r2_h_S1T0 > r2_h_T0S1_min & subset$r2_h_S1T0 < r2_h_T0S1_max,]
  
  if (dim(subset)[1]==0){stop ("No valid observations for the specified ranges of values.")}
  
  mode <- function(data) {
    x <- data
    z <- density(x)
    mode_val <- z$x[which.max(z$y)]
    fit <- list(mode_val= mode_val)
  }
  
  cat("\n\nSummary measures for R2_H (in the subgroup of results where the counterfactual\n")
  cat("correlations fall within prespecified ranges\n")
  cat("##############################################################################\n")
  
  cat("\n\n# R2_H results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")  
  cat("Mean (SD) R2_H: ", format(round(mean(subset$R2_H), 4), nsmall = 4), " (", format(round(sd(subset$R2_H), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(subset$R2_H), 4), nsmall = 4), "; max: ",  format(round(max(subset$R2_H), 4), nsmall = 4), "]", sep="")
  cat("\nMode R2_H: ", format(round(mode(subset$R2_H)$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the R2_H distribution: \n\n")
  quant <- quantile(subset$R2_H, probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)
  
  cat("\n\nNote. The figure is based on ", dim(subset)[1], " observations. \n", sep="")
  
  
  if (Type=="Density"){
    plot(density(subset$R2_H, na.rm = T), xlab=expression(paste(R[H]^2)), ylab="Density", main=MainPlot, lwd=2, ...)
  }
  
  if (Type=="Freq"){
    hist(subset$R2_H, col="grey", main=MainPlot, xlab=expression(paste(R[H]^2)))  
  }
  
  if (Type=="All.Densities"){
    R2_H_No <- subset$R2_H[subset$Monotonicity=="No"]
    R2_H_Surr <- subset$R2_H[subset$Monotonicity=="Surr"]
    R2_H_True <- subset$R2_H[subset$Monotonicity=="True"]
    R2_H_SurrTrue <- subset$R2_H[subset$Monotonicity=="SurrTrue"]
    
    plot(density(subset$R2_H, na.rm = T), xlab=expression(paste(R[H]^2)), ylab="Density", main=MainPlot, lwd=2, col=0, ...)
    
    try(lines(density((R2_H_No), na.rm = T), lty=1, col=1, lwd=3), silent=TRUE) 
    try(lines(density((R2_H_Surr), na.rm = T), lty=2, col=2, lwd=3), silent=TRUE)
    try(lines(density((R2_H_True), na.rm = T), lty=3, col=3, lwd=3), silent=TRUE)
    try(lines(density((R2_H_SurrTrue), na.rm = T), lty=4, col=4, lwd=2), silent=TRUE)
    
    legend(Cex.Position, lwd=c(3, 3, 3, 3), col=c(1, 2, 3, 4), lty=c(1, 2, 3, 4), cex = Cex.Legend,
           legend=c("No monotonicity", "Monotonicity S", "Monotonicity T", "Monotonicity S and T"))
  }
} 