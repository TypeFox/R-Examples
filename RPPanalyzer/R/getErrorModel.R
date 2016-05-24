# Purpose: Fit the error model and append two new columns with the fitted parameters to the data.frame
getErrorModel <- function(dataexpression, verbose=FALSE) {
  
  ## Create data.frame
  
  dataexpression <- dataexpression[!is.na(dataexpression[,1]),]
  
  emptyCols <- rev(which(dataexpression[1,] == ""))[1]
  emptyRows <- rev(which(dataexpression[,1] == ""))[1]
  
  rowNames <- as.character(unlist(dataexpression[1:(emptyRows+1),emptyCols+1]))
  colNames <- as.character(unlist(dataexpression[emptyRows+1, 1:(emptyCols+1)]))
  
  nGpr <- which(rowNames == "gpr")
  nTarget <- which(rowNames == "target")
  nID <- which(colNames == "ID")
  
  slides <- as.character(unlist(dataexpression[nGpr, -(1:(emptyCols+1))]))
  abs <- as.character(unlist(dataexpression[nTarget, -(1:(emptyCols+1))]))
  IDs <- as.character(unlist(dataexpression[-(1:(emptyRows+1)), nID]))
  
  signals <- NULL
  for(i in (emptyCols+2):(dim(dataexpression)[2])) {
    
    description <- dataexpression[-(1:(emptyRows+1)), 1:(emptyCols+1)]
    signal <- as.numeric(as.character(dataexpression[-(1:(emptyRows+1)), i]))
    slide <- as.character(dataexpression[nGpr, i])
    ab <- as.character(dataexpression[nTarget, i])
    
    all <- cbind(slide, ab, description, signal)
    colnames(all) <- c("slide", "ab", colNames, "signal")
    
    signals <- rbind(signals, all)
    
  }
  
  
  ## Subsetting
  
  result <- c()
  ID<-NULL  
  if(verbose) pdf("errorModel.pdf")
  
  for(A0 in unique(signals$slide)) {
    subset1 <- subset(signals, slide==A0)
    
    for(A1 in unique(subset1$ab)) {
      subset2 <- subset(subset1, ab==A1)
      
      
      variances <- t(sapply(as.list(unique(subset2$ID)), function(IDvalue) {
        signals <- subset(subset2, ID==IDvalue)$signal
        variance <- var(signals)
        signal <- mean(signals)
        return(c(signal, variance))}))
      
      variances <- as.data.frame(variances)
      if(verbose) print(variances)
      # - log likelihood function
      ll <- function(var0, varR) {
        y <- variances$V1
        v <- variances$V2
        value <- sum(log(var0 + varR * y^2) + v/(var0+varR*y^2))
        return(value)
      }
      
      gr0 <- function(pars) {
        var0 <- pars
        varR <- 0
        y <- variances$V1
        v <- variances$V2
        DllDvar0 <- sum(1/(var0+varR*y^2) - v/((var0+varR*y^2)^2))
        return(DllDvar0)
      }
      
      grR <- function(pars) {
        var0 <- 0 
        varR <- pars
        y <- variances$V1
        v <- variances$V2
        DllDvarR <- sum((y^2)/(var0+varR*y^2) - v*(y^2)/((var0+varR*y^2)^2))
        return(DllDvarR)
      }
      
      gr <- function(pars) {
        var0 <- pars[1]
        varR <- pars[2]
        y <- variances$V1
        v <- variances$V2
        DllDvar0 <- sum(1/(var0+varR*y^2) - v/((var0+varR*y^2)^2))
        DllDvarR <- sum((y^2)/(var0+varR*y^2) - v*(y^2)/((var0+varR*y^2)^2))
        return(c(DllDvar0, DllDvarR))
      }
      
      
      var0 <- mean(variances$V2[order(variances$V1)[1:(dim(variances)[1]/20)]])
      varR <- 0.03^2 
      
      
      
      fit  <- try(mle(ll, start=list(var0 = var0, varR = varR), gr=gr), silent=TRUE)
      if(class(fit)=="try-error") {
        fit <- try(mle(ll, start=list(var0 = var0, varR = varR), fixed=list(varR=0), gr=gr0), silent=TRUE) 
      } else if(class(fit)!="try-error" & coef(fit)[1] < 0) {
        fit <- try(mle(ll, start=list(var0 = var0, varR = varR), fixed=list(var0=0), gr=grR), silent=TRUE)
      }
      
      if(verbose) print(fit)
      
      
      if(verbose) plot(variances[,1], variances[,2], xlab="Signal", ylab="Variance", log="y")
      y <- seq(min(variances$V1), max(variances$V1), len=1000)
      if(class(fit)!="try-error") {
        var0 <- coef(fit)[1]
        varR <- coef(fit)[2]
        v <- var0 + varR*y^2
        if(verbose) matplot(y, v, type="l", col=2, add=TRUE)
        if(verbose) title(paste(A0,A1))
      }
      
      
      if(class(fit)!="try-error") result <- rbind(result, cbind(subset2, data.frame(var0=coef(fit)[1], varR= coef(fit)[2]))) else result <- rbind(result, cbind(subset2, data.frame(var0=NA, varR= NA)))
      
      
    }
    
    
    
    
    
  }
  
  
  if(verbose) dev.off()
  return(result)
  
  
  
}
