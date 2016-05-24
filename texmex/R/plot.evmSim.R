plot.evmSim <-
function(x, which.plots=1:3, density.adjust=2, print.seed = FALSE , ...){

    varnames <- names(coef(x))

    if (1 %in% which.plots){
        kdes <- list()

        for (i in 1:length(varnames)){
            kdes[[i]] <- density(x$param[, i], n=200, adjust=density.adjust)
        }

        for(i in 1:length(kdes)){
            plot(kdes[[i]], type = "l" , xlab = varnames[i] , ylab="Density" ,
                 main = paste("Posterior for", varnames[i]))
                 polygon(c(rev(kdes[[i]]$x), kdes[[i]]$x),
                         c(rep(0, length(kdes[[i]]$y)) , kdes[[i]]$y),
                         col=6)
                m <- mean(x$param[, i])
    	        ym <- kdes[[i]]$y[abs(kdes[[i]]$x - m) == min(abs(kdes[[i]]$x - m))]
                segments(x0=m, y0=0, y1=ym , x1=m)
    		if (print.seed)
    			title(sub = paste(c("Seed: ", x$seed) , collapse = " "),
                      adj=0)
        }
    } # Close if (1 %in% which plots

	if( 2 %in% which.plots ){
		# Trace plots and moving averages
		for (i in 1:length(varnames)){
    		plot(x$chains[ , i], type = "l" , xlab = "Step number", ylab = paste(varnames[i], "& cumulative mean"))
    		lines(cumsum(x$chains[ , i]) / (1:nrow(x$chains)), lwd=3, col=8)
    		axis(4)
    		abline(v = x$burn, lty=4, lwd=3, col=2)
    		if (print.seed)
    			title(sub = paste(c("Seed: ", x$seed) , collapse = " "), adj=0)
        }
	}

	if(3 %in% which.plots){
        # Plot ACFs

        for(i in 1:length(varnames)){
    		acf(x$param[ , i] , main = paste("ACF for", varnames[i], "\n(thin =",x$thin,")"))
    		if (print.seed)
    			title(sub = paste(c("Seed: ", x$seed) , collapse = " "), adj=0)
        }
    }
	invisible()
}

