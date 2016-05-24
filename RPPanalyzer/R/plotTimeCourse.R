`plotTimeCourse` <- function (x
    ,tc.identifier=c("sample","stimulation","inhibition","stim_concentration")
    ,tc.reference=NULL
    ,plot.split="experiment",file=NULL
    ,arrays2rm=c("protein","Blank"),plotformat="stderr"
    ,log=TRUE,
    color=NULL, xlim=NULL, ylim=NULL) {

    ## select the measurements in the data
    data <- select.measurements(x)

    ## remove arrays without (biological) information
    data <- remove.arrays(data,param="target",arrays2rm=arrays2rm)

    ## create a column to identify the time courses that should printed in one plot
    datax <- create.ID.col(data,sample.id=tc.identifier)

    ### if a reference is given
    #if(!is.null(y)) {

    ### select measurements and remove the control arrays
    #reference <- select.measurements(y)
    #reference <- remove.arrays(reference, param="target", arrays2rm=arrays2rm)

    ## check if the same arrays are available in the reference
    #m <- match(colnames(x[[1]]), colnames(reference[[1]]), nomatch=0)

    #if (any(m==0)) {
    #stop("Some of the arrays in the data are not available in the reference!")
    #} 

    #}

    plotcount <- unique(datax[[4]][,plot.split])
    
#~     yt <- (!is.null(ylim) && ylim=="fixed")
#~     xt <- (!is.null(xlim) && xlim=="fixed")
#~     if(yt || xt) {
#~         lims <- get_ranges(plotcount, plot.split, datax, tc.reference, log)
#~     }
    if(!is.null(file)) {
		pdf(file=file)
		
	}

    for ( i in seq(along=plotcount)){
#~         if(yt || xt) {
#~             if(yt)
#~                 ylim <- lims$ylims[[i]]
#~             if(xt)
#~                 xlim <- lims$xlims[[i]]            
#~         }
        
        #tempdat <- select.sample.group(datax,param=plot.split,sel=plotcount[i])
        groupFactor <- list(plotcount[i])
        names(groupFactor) <- plot.split
        tempdat <- select.sample.group(datax, params=groupFactor)


        ## order data subset
        order.tc <- order(tempdat[[4]][,"identifier"],tempdat[[4]][,"time"])
        tempdat[[1]] <- tempdat[[1]][order.tc,]
        tempdat[[2]] <- tempdat[[2]][order.tc,]
        tempdat[[4]] <- tempdat[[4]][order.tc,]
        


        # create the identifiers for the different time series in one plot

        time.series <- unique(tempdat[[4]][,"identifier"])
    
        # create a nice looking color vector
        if(is.null(color)) {
            color <- rainbow(length(time.series))
        }

        # iterate over all arrays/proteins
        for (k in 1:ncol(tempdat[[1]])){



            # if a reference is given
            if(!is.null(tc.reference)) {

                # get the ID for the reference time series
                refID <- paste(tc.reference, collapse="")

                # search the time series in the data
                refLines <-  tempdat[[4]][,"identifier"]==refID

                # get the zero time point of the time series
                refLines & tempdat[[4]][,"time"]==0

                # check there data available, that means if the combination of parameters was valid for this split
                # if it was, than calculate the reference value
                # otherwise set the reference value to zero or one, depending on the log status
                if(!any(refLines) && log) {

                    ref <- 0

                }
                else if(!any(refLines) && !log) {

                    ref <- 1

                }
                else {

                    ref <- median(tempdat[[1]][refLines & tempdat[[4]][,"time"]==0,k], na.rm=TRUE)

                }
            }

	    # divide by the reference
	    if(!is.null(tc.reference) && log) {

                    tempdat[[1]][, k] <- tempdat[[1]][, k] - ref
                    tempdat[[2]][, k] <- tempdat[[2]][, k] - ref



	    }
	    else if(!is.null(tc.reference) && !log) {

                    tempdat[[1]][, k] <- tempdat[[1]][, k]/ref
                    tempdat[[2]][, k] <- tempdat[[2]][, k]/ref

	    }

            par(lwd=2)

            if(is.null(ylim)) {
                tdk_sd <- sd(tempdat[[1]][,k])
                ylim2 <- c(min(tempdat[[1]][,k])-tdk_sd,max(tempdat[[1]][,k])+tdk_sd)
            } else {
                ylim2 <- ylim
            }
            if(is.null(xlim)) {
                xlim2 <- c(0,max(tempdat[[4]][,"time"]))
            } else {
                xlim2 <- xlim
            }
            # create the plotting area with a simple plot
            plot(tempdat[[4]][,"time"],tempdat[[1]][,k]
                    ,type="n",main=c(paste("time course: ",plotcount[i]),tempdat[[3]]["target",k])
                    ,ylab="signal",xlab="time"
                    ,ylim=ylim2
                    ,xlim=xlim2)
            
            for (j in seq(along=time.series)) {


                time.lines <- which(tempdat[[4]][,"identifier"]==time.series[j])


                expr.data <- cbind(tempdat[[4]][time.lines,"time"]
                        ,tempdat[[1]][time.lines,k]
                        ,tempdat[[2]][time.lines,k]
                        )


                colnames(expr.data) <- c("time",tempdat[[3]]["target",k],"error")
                expr.data <- expr.data[order(expr.data[,"time"]),]


                #if (plotformat=="rawdata" | plotformat=="both" | plotformat=="errbar"){
                if(plotformat %in% c("rawdata", "both", "errbar")) {
#browser()
					if(plotformat=="errbar") {
						errbar.col <- color[j]
					} else {
						errbar.col <- "black"
					}
                    errbar(expr.data[,"time"],expr.data[,2]
                            ,expr.data[,2]+expr.data[,3]
                            ,expr.data[,2]-expr.data[,3]
                            ,col=color[j],add=TRUE,
                            errbar.col=errbar.col)
					if(!plotformat %in% "errbar") {
						lines(expr.data[,"time"],expr.data[,2],col=color[j],lty="dashed",lwd=1)
					}
                }

                #if (plotformat=="spline" | plotformat=="both" | plotformat=="spline_noconf" | plotformat=="errbar"){
                if(plotformat %in% c("spline", "both", "spline_noconf", "errbar", "stderr")) {
                    y <- expr.data[,2]
                    tp <- expr.data[,1]
                    splinemodel <- gam(y~s(tp))
                    # x-vector for the predictions
                    xn=seq(0,max(tp),0.1)
                    ## standard error prediction
                    pred2 <- predict(splinemodel, se.fit=TRUE)
                    if(!plotformat %in% c("spline_noconf", "errbar", "stderr")) {
                        ## generate confidence bands
                        upper <- pred2$fit + pred2$se.fit
                        lower <- pred2$fit - pred2$se.fit
                        band.color <- paste(col2hex(color[j]),"22", sep="")
                        polygon(c(tp, rev(tp)), c(upper, rev(lower)), col=band.color, border=NA)
                        lines(tp, upper, col=paste(col2hex(color[j]),"88", sep=""), lty=2)
                        lines(tp, lower, col=paste(col2hex(color[j]),"88", sep=""), lty=2)
                    }
                    ## add spline to plot
                    pred <- predict(splinemodel,newdata=data.frame(tp=xn))
                    lines(xn, pred, col=color[j], lty="solid", lwd=2)
                    
                    if(plotformat %in% "stderr") {
						 errbar(tp,pred2$fit
                            ,pred2$fit+pred2$se.fit
                            ,pred2$fit-pred2$se.fit
                            ,col=color[j],add=TRUE,
                            errbar.col=color[j])
					}
                    
                   
                }
            }

            legend("topleft",col=color,legend=time.series
                    ,pch=rep(1,length(time.series))
                    ,cex=0.5,lwd=1)

        }
    }
    if(!is.null(file)) {
		dev.off()
	}
}
#~ 
#~ get_ranges <- function(plotcount, plot.split, datax, tc.reference, log) {
#~     ## if ylim="fixed" or ylim=="fixed" get the ranges of x and y axis
#~     ylims <- xlims <- list()
#~    for ( i in seq(along=plotcount)){
#~         #tempdat <- select.sample.group(datax,param=plot.split,sel=plotcount[i])
#~         groupFactor <- list(plotcount[i])
#~         names(groupFactor) <- plot.split
#~         tempdat <- select.sample.group(datax, params=groupFactor)
#~         ## order data subset
#~         order.tc <- order(tempdat[[4]][,"identifier"],tempdat[[4]][,"time"])
#~         tempdat[[1]] <- tempdat[[1]][order.tc,]
#~         tempdat[[2]] <- tempdat[[2]][order.tc,]
#~         tempdat[[4]] <- tempdat[[4]][order.tc,]
#~         # create the identifiers for the different time series in one plot
#~         time.series <- unique(tempdat[[4]][,"identifier"])
#~         # iterate over all arrays/proteins
#~         for (k in 1:ncol(tempdat[[1]])){
#~             # if a reference is given
#~             if(!is.null(tc.reference)) {
#~                 # get the ID for the reference time series
#~                 refID <- paste(tc.reference, collapse="")
#~                 # search the time series in the data
#~                 refLines <-  tempdat[[4]][,"identifier"]==refID
#~                 # get the zero time point of the time series
#~                 refLines & tempdat[[4]][,"time"]==0
#~                 # check there data available, that means if the combination of parameters was valid for this split
#~                 # if it was, than calculate the reference value
#~                 # otherwise set the reference value to zero or one, depending on the log status
#~                 if(!any(refLines) && log) {
#~                     ref <- 0
#~                 }
#~                 else if(!any(refLines) && !log) {
#~                     ref <- 1
#~                 }
#~                 else {
#~                     ref <- median(tempdat[[1]][refLines & tempdat[[4]][,"time"]==0,k], na.rm=TRUE)
#~                 }
#~             }
#~             # divide by the reference
#~             if(!is.null(tc.reference) && log) {
#~                         tempdat[[1]][, k] <- tempdat[[1]][, k] - ref
#~                         tempdat[[2]][, k] <- tempdat[[2]][, k] - ref
#~             }
#~             else if(!is.null(tc.reference) && !log) {
#~                         tempdat[[1]][, k] <- tempdat[[1]][, k]/ref
#~                         tempdat[[2]][, k] <- tempdat[[2]][, k]/ref
#~             }
#~         }
#~         # special treatment of fixed argument for ylim
#~         # use total minimum and maximum as y-axis limits
#~         ylim <- range(tempdat[[1]])
#~         xlim <- range(as.numeric(unique(tempdat[[4]][,"time"])))
#~         ylims[[i]] <- ylim
#~         xlims[[i]] <- xlim
#~     }
#~     return(list(xlims=xlims, ylims=ylims))
#~ }
