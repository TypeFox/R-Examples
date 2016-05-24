
## print, summary and plot method for objects of the class "termstrc_cs"

print.termstrc_cs <- function(x,...) {
  cat("---------------------------------------------------\n")
  if(x$rse) cat("Estimated parameters and robust standard errors:\n") else
  cat("Estimated parameters and standard errors:\n") 
  cat("---------------------------------------------------\n")
  for(i in seq(x$n_group)) {
    print(paste(names(x$alpha)[[i]],":",sep=""))
    cs_coef <- if(x$rse) coeftest(x$regout[[i]],varcovvar=vcovHAC.default) else summary(x$regout[[i]])
    if(x$rse) rownames(cs_coef) <- paste("alpha",c(seq_along(x$alpha[[i]]))) else
     rownames(cs_coef$coefficients) <- paste("alpha",c(seq_along(x$alpha[[i]])))
   print(cs_coef)
   cat("\n")
  x
  }
 }
 

summary.termstrc_cs <-
    function(object,...) {
    x <- object
    RMSE_p <- mapply(function(i) rmse(x$p[[i]],x$phat[[i]]),seq(x$n_group))
    AABSE_p <- mapply(function(i) aabse(x$p[[i]],x$phat[[i]]),seq(x$n_group))
    RMSE_y <- mapply(function(i) rmse(x$y[[i]][,2]*100,x$yhat[[i]][,2]*100),seq(x$n_group))
    AABSE_y <- mapply(function(i) aabse(x$y[[i]][,2]*100,x$yhat[[i]][,2]*100),seq(x$n_group))
    gof <- rbind(RMSE_p,AABSE_p,RMSE_y,AABSE_y)
    colnames(gof) <- names(x$p)
    rownames(gof) <- c("RMSE-Prices","AABSE-Prices","RMSE-Yields (in %)","AABSE-Yields (in %)")
    sumry <- list(gof)
    names(sumry) <- c("gof")
    class(sumry) <- "summary.termstrc_cs"
    sumry
} 


print.summary.termstrc_cs <-
    function(x,...) {
    cat("---------------------------------------------------\n")
    cat("Goodness of fit:\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    print.default(format(x$gof,digits=6,scientific=FALSE),quote=FALSE)
    cat("\n")
    x$gof
    
}


plot.termstrc_cs <-
  function(x,matrange =c(min(mapply(function(i) min(x$y[[i]][,1]), seq(x$n_group))),
                        max(mapply(function(i) max(x$y[[i]][,1]), seq(x$n_group)))),
                        multiple=FALSE, ctype="spot",
                        lwd=2,lty=1,type="l",errors="none",inset=c(0.1,0.3),ask=TRUE, ...) {
    
   # min and max maturity of all bonds in the sample 
     samplemat <- c(min(mapply(function(i) min(x$y[[i]][,1]), seq(x$n_group))),
                    max(mapply(function(i) max(x$y[[i]][,1]), seq(x$n_group))))     
 
    
    cdata <- switch(ctype, "spot" = x$spot,
    					   "forward" = x$forward,
    					   "discount" = x$discount
    					    )
    					   
    cname <- switch(ctype, "spot" = "Zero-coupon yield curve",
    					   "forward" = "Forward rate curve",
    					   "discount" = "Discount factor curve" )
    
    
    # plot all interst rate curves together
    if (multiple) {
    
    plot(x=cdata,multiple=multiple, expoints=NULL,lwd=lwd,type=type,...) }
  
	 if (!multiple && ctype %in% c("spot", "forward", "discount")){
        old.par <- par(no.readonly = TRUE)
        if(x$n_group !=1) par(ask=ask)
        
    	# plot each interest rate curve seperately
    	for (k in seq(x$n_group)  ) 
    	{
    	
    	plot.ir_curve(cdata[[k]], ylim=c(0, max(cdata[[k]][,2]) + 0.01 )*100,
    	xlim=c(max(floor(min(x$y[[k]][,1])),matrange[1]),
             min(ceiling(max(x$y[[k]][,1])),matrange[2])), lwd=lwd,type=type,...)
    	
    	
    	 
    	title(x$group[k])
    	 
    	if(ctype=="spot") {points(x$y[[k]][,1],x$y[[k]][,2]*100,col="red") 
    	  # lower ci         
          lines(cdata[[k]][,1],cdata[[k]][,3]*100, type="l", lty=3, col="steelblue" )   
          # upper ci 
          lines(cdata[[k]][,1],cdata[[k]][,4]*100, type="l", lty=3, col="steelblue")
    	  # knot points 
    	  abline(v=c(x$knotpoints[[k]]),lty=2, col="darkgrey")
            legend("bottom",legend=c("Zero-coupon yield curve",
    	  if(x$rse) "95 % Confidence interval (robust s.e.)" else "95 % Confidence interval" ,"Yield-to-maturity", "Knot points"),
    	  col=c("steelblue","steelblue","red", "darkgrey"),
    	  lty = c(1,3,-1,2), pch=c(-1,-1,21,-1))
	
    	     } else  legend("bottom",legend=cname,col=c("steelblue"), lty = lty , pch=(-1))

    	
    	 }
        on.exit(par(old.par))
 	}
    	
    # plot spread curves 
    if(ctype == "spread") {plot(x$spread,expoints=NULL,
    	xlim= c(max(floor(samplemat[1]),matrange[1]),
  	     min(ceiling(samplemat[2]),matrange[2],max(mapply(function(i) 
	     max(x$spread[[i]][,1]),seq(x$spread))))),lwd=lwd ,...) 
       }
    						
    						
     # plot errors 
    if(errors %in% c("price", "yield")){
    	
    	edata <- switch(errors,"price" = x$perrors, "yield"= x$yerrors )

        if(x$n_group == 1) ask= FALSE

        for(k in seq(x$n_group)){
     		plot.error(edata[[k]],ask=ask
                ,main=x$group[k],ylab=paste("Error ",paste(errors,"s)",sep=""),sep=" ("),...)
    		
    		legend("bottomright", legend=c(paste("  RMSE",
    		switch(errors,"price" = round(rmse(x$p[[k]],x$phat[[k]]),4),
                       "yield" = round(rmse(x$y[[k]][,2],x$yhat[[k]][,2]),4)) ,sep=": "),
                        paste("AABSE",switch(errors,"price" = round(aabse(x$p[[k]],x$phat[[k]]),4),
                        "yield" = round(aabse(x$y[[k]][,2],x$yhat[[k]][,2]),4)),sep=": ")),bty="n", inset=inset) 
    		
    	  }
    	
     }						
    					     					        							
} 
