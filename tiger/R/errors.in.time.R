errors.in.time <- function(xval, result, solution, rain.data=NULL, show.months=FALSE,
new.order=1:solution, x.range=1:length(xval),
pmax=max(c(result$measured, result$modelled), na.rm=TRUE),
data.colors=data.frame(measured=c("grey"), modelled=c("black"), rain=c("black")),
clusterPalette= rainbow(solution),
color.cut.off=0, frac.max=0.7, frac.min=0.4, grid.nx=0,
legend.pos="topleft", show.data=TRUE, show.errors=TRUE, show.data.model=show.data, show.data.measured=show.data, ...
 ){


    pmax.orig <- pmax
    pmin<-frac.min*pmax
    pmax<-frac.max*pmax
    pdist <- (pmax-pmin)/solution
    old.pal <- palette(clusterPalette)

    if(show.data.model | show.data.measured){
        plot(range(xval[x.range]), range(0,pmax.orig), ylab="discharge/mm/h", type="n", ...)
    }  else {
        plot(range(xval[x.range]), range(0,pmax.orig), ylab="", type="n", yaxt="n", ...)
    }

    if(show.months){
        min.month <- strftime(min(as.POSIXlt(xval[x.range])), "%m.%y")
        my.min <- strptime(paste(1,min.month,sep="."), "%d.%m.%y")
        max.month <- strftime(max(as.POSIXlt(xval[x.range])), "%m.%y")
        my.max <- strptime(paste(1,max.month,sep="."), "%d.%m.%y")
        month.ticks <-seq(my.min,my.max,by="month") 
        month.lab <- strftime(as.POSIXlt(month.ticks), "%b")
        month.lab[strftime(as.POSIXlt(month.ticks), "%m")=="01"] <- NA
        axis(1, at=month.ticks, labels=month.lab, cex.axis=0.4, tcl=-0.15, mgp=c(3,0.2,0))
    }

    leg.col=c()
    leg.text=c()

    if(show.data.model){
        if(result$multi.model){
            for(i in 1:result$count.model){
                reduced.lines(xval[x.range], result$model[i,x.range], type="l", lwd=3, col=as.character(data.colors$modelled))
            }
        } else {
            reduced.lines(xval[x.range], result$model[x.range], type="l", lwd=3, col=as.character(data.colors$modelled))
        }
        leg.col = as.character(data.colors$modelled)
        leg.text = "modelled"
    }
    if(show.data.measured){
        reduced.lines(xval[x.range], result$measured[x.range], type="l", lwd=2.5, col=as.character(data.colors$measured))
        leg.col = c(leg.col, as.character(data.colors$measured))
        leg.text = c(leg.text, "measured")
	if(!is.null(rain.data)){
		old.usr <- par("usr")
		par(usr=c(old.usr[1:2],  max(rain.data)*2, -0.1))
		lines(xval[x.range], rain.data[x.range], type="h")
		axis(side=4)
		mtext("rain/mm/h", line=2, side=4, cex=0.5)
		par(usr=old.usr)
	}
    }
    if(show.data.model | show.data.measured){
        legend(legend.pos, inset=0.05, lty=c(1,1), lwd=c(3,2.5),col=leg.col,  legend=leg.text)
    }

    if(show.errors){
        region.proportion <- cluster.proportion(result=result, solution=solution, new.order=new.order)
        #redefine x.range according to step size
	x.range.orig <- x.range
	if(class(x.range)=="logical") x.range <- which(x.range)
        x.range <- seq(min(x.range), max(x.range), by=result$step.size)
	if(NCOL(region.proportion)>length(x.range)) region.proportion <- region.proportion[,x.range]
        stopifnot(length(x.range)==dim(region.proportion)[2])
        min.diff <- min(diff(xval[x.range]))
        try(units(min.diff) <- "secs", silent=TRUE)
        min.diff <- as.numeric(min.diff)
        myFrom <- xval - min.diff/2
        myTo <- xval + min.diff/2

        for(the.region in 1:solution){
           currmin <- pmax-(the.region-1)*pdist
           currmax <- pmax-the.region*pdist
           cols <- color.factor(color=palette()[the.region],
                                value= region.proportion[the.region,],
                                max=1)
           cols[region.proportion[the.region,]<=color.cut.off] <- NA
           rect(myFrom[x.range],currmin, myTo[x.range], currmax, col=cols, border="transparent")

           #lines(xval[sel][subsel], result$measured[sel][subsel], type="l", col=the.region)
        }
        axis(side=4, at= seq((pmax-pdist/2), pmin+pdist/2, length.out = solution), labels=LETTERS[1:solution], mgp=c(3,0.7,0), las=2, cex.axis=0.6 )
    }
    grid(nx=grid.nx, ny=0)

     palette( old.pal )

}
