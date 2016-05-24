plot.gsym.point <-
function(x, legend=TRUE,...) {
	op <- par(pty = "s")
      methods <- x[x$methods]
	n.levels.cat <- if(is.null(x$levels.cat)) {1} else {length(x$levels.cat)}
	levels.cat <- if(is.null(x$levels.cat)) {NULL} else {x$levels.cat}
      n.plots=0
	for (i in 1:n.levels.cat) 
	{
		for(j in 1:length(methods)) 
		{
			if(length(methods[[j]][[i]][["optimal.result"]][[1]])== 0) 
			{
				if(is.null(x$levels.cat)) {
					cat(paste(names(methods)[j], ": There are no cutoff values that fulfill the criterion \n", sep = ""))
				} else {
					cat(paste(names(methods)[j], ": There are no cutoff values that fulfill the criterion for ", levels.cat[i], "\n", sep = ""))
				}
			} 			
			
     	     		m <- methods[[j]][[i]]			          
                      
                                  
    			if(n.plots > 0)
            	{
     	    			readline("Press return for next page....")
      		}  

                  if (n.levels.cat >1)
			{
				# Marker in the healthy population:
      			X0 <- x$data[x$data[,x$call$status] == x$call$tag.healthy & x$data[,x$call$categorical.cov] == levels.cat[i], x$call$marker]
      			# Marker in the diseased population:
      			X1 <- x$data[x$data[,x$call$status] != x$call$tag.healthy & x$data[,x$call$categorical.cov] == levels.cat[i], x$call$marker]
			}

			else
			{
				# Marker in the healthy population:
      			X0 <- x$data[x$data[,x$call$status] == x$call$tag.healthy, x$call$marker]
      			# Marker in the diseased population:
      			X1 <- x$data[x$data[,x$call$status] != x$call$tag.healthy, x$call$marker]
			}

				
      		x0 = relative_sample(X0,X1,0)
      		rho = m[["rho"]]
                        
			main <- paste("Empirical ROC Curve and line y = 1-",rho, "x", "\n", " ", ifelse(is.null(levels.cat), "", levels.cat[i]), sep = "")
     	     		
      		parameters = c(rho,x0)
      		x_axis = seq(0,1,length.out=1000)
      		y_axis_1 = 1-cdf_empirical_dist(x0,1-x_axis)
      		plot(x_axis,y_axis_1,type="l",xlim=c(0,1),ylim=c(0,1),xlab=expression("1-Sp(c)"),ylab=expression("Se(c)"), main = main, cex.main=1,...)
      		lines(x_axis, 1-rho*x_axis, lty = 2)
                        
      		## Legend for the AUC:
      		if (legend)
      		{
      			legend.text <- paste("AUC: ",round(m[["AUC"]][[1]], 3))
      			legend("bottom", legend.text, bty = "n")	

      		}
                  n.plots=n.plots + 1

		}
	}
	par(op)
}
