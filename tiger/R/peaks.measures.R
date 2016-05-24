peaks.measures <- function(result, show.measures=1:num.measures,
synthetic.peaks.col = c(2:8,2:8),mfrow=c(2,3),col.best.match="black", do.out=rep(TRUE, length(show.measures)), single.errors=FALSE, show.legend=TRUE ,
show.main=TRUE, y.range=NULL){
    #ToDo: Make this method working either with tiger results or results from the synthetic.errors function
    n.level <- dim(result$measures.synthetic.peaks)[2]
    n.errors <- dim(result$measures.synthetic.peaks)[1]
    num.measures <- length(result$measures)

    my.min <- function(x, na.rm, do.out=TRUE){
         return(min(x[!x %in% boxplot.stats(x, do.out=do.out)$out], na.rm=na.rm))
    }
    my.max <- function(x, na.rm, do.out=TRUE){
         return(max(x[!x %in% boxplot.stats(x, do.out=do.out)$out], na.rm=na.rm))
    }
    palette("default")

    top.mar <- ifelse(show.main, 3, 1)
par(mfrow=mfrow,mar=c(2,2,top.mar,1)+0.1, oma=c(0,0,0,0)+0.1)
for(measure in show.measures){
	if(measure==show.measures[length(show.measures)] & single.errors){
		xaxt="s"
		bottom.mar=2
	} else {
		xaxt="n"
		bottom.mar=0
	}
        data <- result$measures.synthetic.peaks[,,measure]
	if(show.main){
		lab <- result$names[[measure]]
		lab <- substitute(" "*a, list(a=lab))
	} else {
		lab <- ""
	}
        best <- result$best.value.location$all.values[measure]
	prepare.plot <- function(xaxt="s", yaxt="s", adjust.range=FALSE){
		if(is.null(y.range)){
			y.range <- c(my.min(c(data,best), na.rm=TRUE, do.out=do.out[which(show.measures %in% measure)]), my.max(c(data,best), na.rm=TRUE, do.out=do.out[which(show.measures %in% measure)]))
			if(adjust.range){
				data.range <- range(data[error,])
				if( data.range[2] < y.range[1] | 
 				    data.range[1] > y.range[2]){
			    		y.range <- data.range
					yaxt="s"
					the.mar <- par()$mar
					the.mar[2] <- 2.1
					par(mar=the.mar)
		    		}
			}
		}
		plot(range(1:n.level), y.range , type="n", main=lab, xlab="", ylab=lab, xaxt=xaxt, yaxt=yaxt)
		lines(c(0,n.level+1),c(best, best) , lwd=2, col=col.best.match)
	}
	if(!single.errors) prepare.plot()
        for(error in 1:n.errors){
		if(single.errors){
			if(error==1){
				yaxt="s"
				par(mar=c(bottom.mar,2,top.mar,1)+0.1)
			}else {
				yaxt="n"
				par(mar=c(bottom.mar,0,top.mar,1)+0.1)
			}
		       	prepare.plot(xaxt=xaxt, yaxt=yaxt, adjust.range=TRUE)
		}
        	points(data[error,], pch=error, col=synthetic.peaks.col[error])
        }
}
	if(show.legend){
		plot(c(0,1),c(0,1), type="n", xaxt="n", yaxt="n", ylab="", bty="n")
		legend("top",result$error.names, pch=1:n.errors, col=synthetic.peaks.col, cex=0.8)
	}
}
