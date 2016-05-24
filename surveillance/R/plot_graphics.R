###################################################
### chunk number 1: 
###################################################

plot.survRes <- function(x, method=x$control$name, disease=x$control$data, startyear = 2001, firstweek = 1, legend=TRUE,...){

        survResObj <- x
        observed <- survResObj$disProgObj$observed[survResObj$control$range]
        state    <- survResObj$disProgObj$state[survResObj$control$range]
        
        # width of the column
        tab <- 0.5

        # left/right help for constructing the columns
        observedxl <- (1:length(observed))-tab
        observedxr <- (1:length(observed))+tab
        upperboundx <- (1:length(survResObj$upperbound))-0.5

        ######################################################################
        # Generate the axis labelling
        ######################################################################
        #Compute how much has been cut off by the specification of range
        start <- min(survResObj$control$range)-1
        startyear <- startyear + start %/% 52
        firstweek <- firstweek + start %%  52

        # get the number of quarters lying in range 
        myat.week <- seq(ceiling((52-firstweek+1)/13) * 13 + 1, length(observed)+(floor((52-firstweek + 1)/13) * 13 +1), by=13)
        # get the right year order
        year <- (myat.week - 52) %/% 52 + startyear
        # function to define the quarter order
        quarterFunc <- function(i) { switch(i+1,"I","II","III","IV")}
        # get the right number and order of quarter labels
        quarter <- sapply( (myat.week-1) %/% 13 %% 4, quarterFunc)
        # get the positions for the axis labels
        myat.week <- myat.week - (52 - firstweek + 1)


        #Find out how much the axes are scaled
        cex <- par()$cex.axis

        # construct the compute axis label
        if (cex == 1) {
          mylabels.week <- paste(year,"\n\n",quarter,sep="")
        } else {
          mylabels.week <- paste(year,"\n",quarter,sep="")
        }

        #Remove NaNs from upperbound (real problem: where do they come from)
        #survResObj$upperbound[is.nan(survResObj$upperbound)] <- 0

        # control where the highest value is
        max <- max(max(observed),max(survResObj$upperbound))


        #Generate the matrices to plot
        xstuff <- cbind(observedxl, observedxr, upperboundx)
        ystuff <-cbind(observed, observed, survResObj$upperbound)
        

        #Plot the results using one Large plot call
        matplot(xstuff,ystuff ,
                t="hhs", lty=c(1,1,1), col=c(1,1,4), ylim=c(-1/20*max, max),
                xlab = "time", ylab = "No. of infected", axes = FALSE)#hoehle, ...)
        if (!is.null(survResObj$aggr)) {
          points(upperboundx+tab,survResObj$aggr,col=1)
        }
        
        for(i in 1:length(observed)){
          matlines( c(i-tab, i+tab), c(observed[i],observed[i]) )
          if(survResObj$alarm[i] == 1)
            matpoints( i, -1/40*max, pch=24, col=2)
          if(state[i] == 1)
            matpoints( i, -1/20*max, pch=24, col=3)
        }
      
        if (disease != "") { disease <- paste("of ",disease," ",sep="") }
        title(paste("Analysis ", as.character(disease), "using ", as.character(method),sep=""))

        # parameters for the legend placement to the right upper side
        xlegpos <- 1/4
        ylegpos <- 1

        # check where to place the legend. If the left upper side is free place it there
        if (max * 2/3 >= max(
                      max(observed[1:floor(1/4 * length(observed))]),
                      max(survResObj$upperbound[1:floor(1/4 * length(survResObj$upperbound))])
                      )) {
          xlegpos <- 0
        }

        if (legend) {
          legend(xlegpos*length(observed)/sqrt(cex), ylegpos*max,
                 legend=c("Infected", "Threshold", "Computed Alarm", "Defined Alarm"),
                 lty=c(1,1,NA,NA), col=c(1,4,2,3), pch=c(NA,NA,24,24),cex=cex)
        }

        axis( at=myat.week , labels=mylabels.week , side=1, line = 1 )
        axis( side=2 )
        
        invisible()
}


