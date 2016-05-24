anm.ls.reg<-function (X, Y, parameter = "slope",  
    nmax=50, interval = 0.1, col = "red",...) 
{
    lin.mod <- lm(Y ~ X)
    b0 <- summary(lin.mod)$coefficients[1]
    b1 <- summary(lin.mod)$coefficients[2]
    stb0 <- summary(lin.mod)$coefficients[3]
    stb1 <- summary(lin.mod)$coefficients[4]
    
    if(parameter == "intercept") {poss <- seq(b0 - (2 * stb0), b0 + (2 * stb0), length.out=nmax)}
    if(parameter == "slope"){poss <- seq(b1 - (2 * stb1), b1 + (2 * stb1), length.out=nmax)}
     
    fitted <- matrix(ncol = length(poss), nrow = length(X))
    if(parameter == "intercept") {for(i in 1:ncol(fitted)) fitted[,i] <- poss[i] + X * b1}  
    if(parameter == "slope"){for(i in 1:ncol(fitted)) fitted[,i] <- b0 + poss[i] * X}      
    
    sq.res <- matrix(nrow = length(X), ncol = length(poss))
    for (i in 1:length(poss)) {
        sq.res[, i] <- (fitted[,i] - Y)^2
    }
        ss.res <- apply(sq.res, 2, sum)
  
    
    
            dev.new(height = 4, width = 8)
            par(mfrow = c(1, 2), mar = c(4.4, 4.5, 1, 0.5), cex = 0.9)
      
               for (i in 1:length(poss)) {
                dev.hold()
                if(i < length(poss)){
                plot(X, Y) 
                abline(b0, b1, lwd = 1.4)
                if(parameter == "intercept") abline(poss[i],b1, col="gray")
                if(parameter == "slope") abline(b0,poss[i], col="gray")
                segments(X, Y, X, fitted[,i],lty = 2, col = col) 
                }
                
                if(i == length(poss)){
                plot(X, Y, xlab = expression(italic(X)), ylab=expression(italic(Y)))
                abline(b0, b1, lwd = 1.4)
                segments(X, Y, X, fitted(lin.mod),lty = 2, col = col)
                }
                          
                plot(poss, ss.res, ylab = "Residual sum of squares", xlab = 
                ifelse(parameter == "intercept","Intercept estimate","Slope estimate"), type = "n") 
                arrows(poss[i], ss.res[i], poss[i + 1], ss.res[i + 1], col = col, length = 0.15, lwd = 1)
                points(poss[1:i], ss.res[1:i], lty = 2, col = col, lwd = 1, type = "l")
                  
                if(i < length(poss)){
                legend("topright", legend = bquote(paste("SS = ", .(ss.res[i]))), cex = 0.9, bty = "n")
                }
                
                 if(i == length(poss)){
                legend("topright", legend = bquote(paste("SS = ", .(sum(resid(lin.mod)^2)))), cex = 0.9, bty = "n")
                }
                
                dev.flush()
                Sys.sleep(interval)
                }
                     
           abline(v=ifelse(parameter=="intercept",b0,b1),lty=2)
           if(parameter=="intercept") legend("center",legend = bquote(paste(hat(beta[0])," = ",.(b0))),box.col="white", bg = "white",)
           if(parameter=="slope") legend("center",legend = bquote(paste(hat(beta[1])," = ",.(b1))),box.col="white", bg = "white")
          
}