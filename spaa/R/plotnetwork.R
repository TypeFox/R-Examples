plotnetwork <-
function(datainput, interval = 8, 
         xlim = c(-2.5, 5), ylim=c(-3.2,3.2), 
		 lty = c(1,2,3,4,4,3,2,1,5),value = "r", 
		 legendx = 3, legendy = 0, right = 1.2,
         intcept = 0.22, left = 0.35, linelength = 0.3, cex = 3,
         lwd = 1.5, show.legend = TRUE, digits = 2, dit = 1.2, 
         number.label = TRUE, text.label = TRUE, 
		 linecol = c("red", "black"), ...){
    datainput <- as.matrix(datainput)
    if (!nrow(datainput)==ncol(datainput)){
        stop("The input data must be a square matrix. \n")
    }
    if(!length(lty)== interval + 1){
	    stop("Number of line styles should be equal to interval +1. \n")
	}
    if(nrow(datainput)>15){
	    cat("Warning: Too many points too draw, try to simplify the input data. \n")
	}
    if (!is.matrix(datainput)){
       if(!is.numeric(datainput)){
		 stop("Only numberic elements can be included in the correlation matrix. \n")
		 }
    }
     if (digits > 2){
	     cat("Warning: Digits of value are suggested to be 1 or 2.")
	 }
    npoints = nrow(datainput)
    if (npoints > 10){
	    cat("Warning: Relationship between more than 10 sites \n
            may produce very complicated conections and difficult to explain. \n")
	}
  
    plot(0,0, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, xlab="", ylab="")
    interval = interval
    datainputd = as.dist(datainput)
    seqn  <- (seq(0, interval^2, by = interval))/(interval^2)
    limit0 <- round((min(datainputd) + (max(datainputd)-min(datainputd))*seqn), digits)
    limit <- limit0[c(-1, -length(limit0))]
    limit <- sort(limit, decreasing = TRUE)
    
    subangle <- 2*pi/npoints
    r=2
    xdep <- c(); ydep <- c(); txdep <- c(); tydep <- c()
    for (k in 1:npoints){ 
        r = r
        x = r*sin(subangle * k)
        y = r*cos(subangle * k)
        tx = dit*x
        ty = dit*y
        xdep <- append(xdep, x)
        ydep <- append(ydep, y)
        txdep <- append(txdep, tx)
        tydep <- append(tydep, ty)
    }
  
    for (i in 1:length(xdep)){
        for (j in 1:length(ydep)){
            if(datainput[i,j] >  0) color = linecol[1]
            if(datainput[i,j] <= 0) color = linecol[2]
            if (datainput[i,j] >= limit[1]){
		        segments(x0 = xdep[i], y0 = ydep[i], x1 = xdep[j], y1 = ydep[j], 
				         lty = lty[1], col= color, lwd = abs(datainput[i,j])*4)}
                    if (datainput[i,j] <  limit[length(limit)]){ 
	                      segments(x0 = xdep[i], y0 = ydep[i], x1 = xdep[j], 
						     y1 = ydep[j], lty = lty[length(limit)+1], 
							 col= color, lwd = abs(datainput[i,j])*4)
						}else{
                              for (n in 2:length(limit)){ 
                                  if (datainput[i,j] <  limit[n-1] & datainput[i,j] >= limit[n]){
                                       if (datainput[i,j] > 0) color = linecol[1]
                                       if (datainput[i,j] < 0) color = linecol[2]
                                       segments(x0 = xdep[i], y0 = ydep[i], x1 = xdep[j], y1 = ydep[j], lty = lty[n], col = color, lwd = abs(datainput[i,j])*4)
                                     }
                                 }
                              }
        
		
		
		
		}
    }
   
   if(number.label){
       points(xdep, ydep, cex=cex, pch = 21, bg = "white", lwd = lwd,)
       text(xdep, ydep, 1:npoints)
   }
   if(text.label){
       text(txdep, tydep, colnames(datainput)) 
   }
   
   if(show.legend){
       for (k in 1:length(limit)){
           if (limit[k] > 0){
	           col = linecol[1]
			   }else{
			           col = linecol[2]
				    }
                   if (k == 1){
                     text (legendx + right, legendy, formatC(paste(value ,">=",
         					 sprintf("%.2f",limit[k])), width =10), )
                     segments(legendx - left, legendy, legendx- left + 
					         linelength, legendy, lty = lty[k], 
							 col = col, lwd = limit[k]*4)
                    }
                   if (k == length(limit)){
                     text (legendx + right, legendy-intcept*k, 
					      formatC(paste(value ,"<", sprintf("%.2f",
						  limit[k])), width = 10))
                     segments(legendx- left, legendy-intcept*k, 
					      legendx- left + linelength, legendy-intcept*k, 
						  lty = lty[k+1], col = col, lwd = abs(limit[k])*4)
                    }
                   else{
                        text (legendx + right, legendy-intcept*k, 
					    formatC(paste(sprintf("%.2f",limit[k+1]),"=<",value ,"<", sprintf("%.2f",limit[k])),width = 10))
                        segments(legendx- left, legendy-intcept*k, legendx- left + linelength, legendy-intcept*k, lty = lty[k+1], col = col, lwd = abs(limit[k])*4)
                      }

	     }
    }
}

