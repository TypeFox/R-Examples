bbeta.output <- function (x, plot = FALSE, ...){
    
    FF <- x$Theta
    cat("\n", "Posterior density estimate:", "\n")
    overall <- rbind(apply(FF, 2, mean), apply(FF, 2, sd))
    rownames(overall) <- c("Post. mean:", "Post. sd:")
    colnames(overall) <- paste("theta", seq(1, x$nstats), " (", 
        x$specs[seq(1, x$nstats)], ")", sep = "")
    all <- as.table(overall)
    print(overall)
    
    cat("\n","Acceptance rate:", x$acc.rate,"\n")
    
    if(plot==TRUE){
         dev.new()
         K <- mcmc(data = FF)
         par(mfrow = c(min(4, x$nstats), 3), 
             oma = c(0, 0, 3, 0), 
             mar = c(4, 3, 1.5, 1))
         for (i in 1:x$nstats) {
             if (i %in% c(5, 9, 13)) {
                 dev.new()
                 par(mfrow = c(min(4, x$nstats - (i - 1)), 3), 
                     oma = c(0, 0, 3, 0), 
                     mar = c(4, 3, 1.5, 1))
             }
             plot(density(FF[, i]), 
                  main = "", 
                  axes = FALSE, 
                  xlab = bquote(paste(theta[.(i)], " (", .(x$specs[i]), ")")), 
                  ylab = "", 
                  lwd = 2)
             axis(1)
             axis(2)
             plot(FF[, i], type = "l", xlab = "Iterations", ylab = "")
             autocorr.plot(K[, i], auto.layout = FALSE, ...)
             if (i %in% union(x$nstats, seq(4, x$nstats, 4)))
                 title(paste("MCMC output for Model: y ~", x$formula[3]), outer = TRUE)  
		 }	
	}
}
