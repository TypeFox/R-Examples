hpbayes.plot <-
function(..., nrisk = NULL, ndeath = NULL, age, hpp, yrange = c(0, 0.8), xrange = c(0, 85), log = FALSE, plotdata = FALSE, plotpost = TRUE, data.type ="b", post.type = "l", line.col=c("grey", "blue", "red", "dark green"), CI=95) 
{
	par(...)
	lc <- as.matrix(line.col)
	loCI <- ((100-CI)/2)/100
	hiCI <- 1-(((100-CI)/2)/100)
	
    hpq <- hp.nqx(H.out = hpp, age = age)
    hpq.med <- rep(NA, length(age))
    for (i in 1:length(hpq.med)) {
        hpq.med[i] <- median(hpq[, i])
    }
    hpq5 <- rep(NA, length(age))
    for (i in 1:length(hpq5)) {
        hpq5[i] <- quantile(hpq[, i], probs=loCI)
    }
    hpq95 <- rep(NA, length(age))
    for (i in 1:length(hpq95)) {
        hpq95[i] <- quantile(hpq[, i], probs=hiCI)
    }
    if (log == FALSE) {
    	if(plotpost==TRUE){
        plot(age, hpq[1, ], type = post.type, lty = 1, lwd = 0.5, 
            col = lc[1,], ylab = expression(paste(""[n], q[x])), 
            xlab = "age", ylim = yrange, xlim = xrange)
        for (i in 2:nrow(hpq)) {
            points(age, hpq[i, ], type = post.type, lty = 1, 
                lwd = 0.5, col = lc[1,])
        }
        
        if (plotdata == TRUE) {
            q <- ndeath/nrisk
            points(age, q, type = data.type, col = lc[4,], 
                cex = 0.6)
        }
        
        points(age, hpq.med, type = post.type, lty = 1, col = lc[2,], 
            lwd = 0.6)
        points(age, hpq5, type = post.type, lty = 2, col = lc[3,], 
            lwd = 0.5)
        points(age, hpq95, type = post.type, lty = 2, col = lc[3,], 
            lwd = 0.5)
    }
		if(plotpost==FALSE){
			if(plotdata==TRUE) {
				q <- ndeath/nrisk
				plot(age, q, type=data.type, col = lc[4,],
					ylab = expression(paste(""[n], q[x])), xlab = "age", ylim = yrange, xlim = xrange)
				points(age, hpq.med, type = post.type, lty = 1, col = lc[2,], lwd = 0.6)
        		points(age, hpq5, type = post.type, lty = 2, col = lc[3,], lwd = 0.5)
       	 		points(age, hpq95, type = post.type, lty = 2, col = lc[3,], lwd = 0.5)
				}
			if(plotdata==FALSE) {
				plot(age, hpq.med, type = post.type, lty = 1, col = lc[2,], lwd = 0.6, 
					ylab = expression(paste(""[n], q[x])), xlab = "age", ylim = yrange, xlim = xrange)
				points(age, hpq5, type = post.type, lty = 2, col = lc[3,], lwd = 0.5)
       	 		points(age, hpq95, type = post.type, lty = 2, col = lc[3,], lwd = 0.5)
				}
    }
    }
    if (log == TRUE) {
        q <- ndeath/nrisk
        lq <- log(q)
        lhpq <- log(hpq)
        lhpq.med <- log(hpq.med)
        lhpq5 <- log(hpq5)
        lhpq95 <- log(hpq95)
        yrange <- c(min(lhpq, lq), max(lhpq, lq))
        xrange <- xrange
    	
    	if(plotpost==TRUE){
        plot(age, lhpq[1, ], type = post.type, lty = 1, lwd = 0.5, 
            col = lc[1,], ylab = expression(log(paste(""[n], q[x]))), 
            xlab = "age", ylim = yrange, xlim = xrange)
        for (i in 2:nrow(lhpq)) {
            points(age, lhpq[i, ], type = post.type, lty = 1, 
                lwd = 0.5, col = lc[1,])
        }
        
        if (plotdata == TRUE) {
            q <- ndeath/nrisk
            points(age, lq, type = data.type, col = lc[4,], 
                cex = 0.6)
        }
        
        points(age, lhpq.med, type = post.type, lty = 1, col = lc[2,], 
            lwd = 0.6)
        points(age, lhpq5, type = post.type, lty = 2, col = lc[3,], 
            lwd = 0.5)
        points(age, lhpq95, type = post.type, lty = 2, col = lc[3,], 
            lwd = 0.5)
    }
		if(plotpost==FALSE){
			if(plotdata==TRUE) {
				plot(age, lq, type=data.type, col = lc[4,],
					ylab = expression(log(paste(""[n], q[x]))), xlab = "age", ylim = yrange, xlim = xrange)
				points(age, lhpq.med, type = post.type, lty = 1, col = lc[2,], lwd = 0.6)
        		points(age, lhpq5, type = post.type, lty = 2, col = lc[3,], lwd = 0.5)
       	 		points(age, lhpq95, type = post.type, lty = 2, col = lc[3,], lwd = 0.5)
				}
			if(plotdata==FALSE) {
				plot(age, lhpq.med, type = post.type, lty = 1, col = lc[2,], lwd = 0.6, 
					ylab = expression(log(paste(""[n], q[x]))), xlab = "age", ylim = yrange, xlim = xrange)
				points(age, lhpq5, type = post.type, lty = 2, col = lc[3,], lwd = 0.5)
       	 		points(age, lhpq95, type = post.type, lty = 2, col = lc[3,], lwd = 0.5)
				}
    }
}
}
