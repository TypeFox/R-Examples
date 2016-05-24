# Cloud Utilization Plot

cloudUtilPlot <- function(group, 
    begin,
    end,
    id,
    colormap=rainbow(length(unique(group))), 
    normalize=TRUE,
    plotConcurrent=TRUE,
    plotConcurrentMax=FALSE,
    main="") {


    # determine the number of compute runs/replicates or whatsoever
    cloudClass<-sort(unique(group))

    if ( normalize == FALSE ) {
        if (nchar(main) == 0){
            main<-"normalized utilization plot"
        }
        plot(as.POSIXct(end,origin="1970-01-01"),
            id,
            main=main,
            type="n",
            xlim=range(c(begin, end), na.rm = T ),
            ylab=list("job id"),
            xlab=list("time [POSIXct since 1970-01-01]"))
    } else {
        # compute the xlim for the plot function
        if (nchar(main) == 0){
            main<-"absolut utilization plot"
        }
        max.g<-rep(0, length(cloudClass))
        for (i in 1:length(cloudClass)){
            minBegin <- min(begin[ group == cloudClass[i] ], na.rm = T)
            s.c.begin <- begin[ group == cloudClass[i] ] 
            s.c.end <- end[ group == cloudClass[i] ] 
            max.g[i] <- max(s.c.end-minBegin, na.rm=T )
        }
        plot(begin,
            id,
            main="", 
            xlim=c(0,max(max.g, na.rm=T)),
            type="n",
            ylab=list("job id"),
            xlab=list("time [in seconds]"))
    }

    # plot cloud utilization plot
    for (i in 1:length(cloudClass)){
        s.c<-group[ group == cloudClass[i] ]

        if ( normalize == T ) {
            minBegin    <- min(begin[ group == cloudClass[i] ], na.rm = T)
        } else {
            minBegin <- 0
        }

        s.c.begin   <- begin[ group == cloudClass[i] ] - minBegin
        s.c.end     <- end[ group == cloudClass[i] ] - minBegin
        s.c.id      <- id [ group == cloudClass[i] ]

        segments(s.c.begin,s.c.id, s.c.end, s.c.id,col=colormap[i])

        abline(v=max(s.c.end), col=colormap[i])
        text(max(s.c.end), 0.6*max(id), cloudClass[i], cex=0.75,srt=90)

        if ( plotConcurrent==T & normalize == T){
            t<-0;
	        n<-100;
	        tt<-rep(0,n)
	        jj<-rep(0,n)
        	idx<-round(seq(0,max(s.c.end, na.rm = T),length=n))
	        for (ii in idx){
		        t<-t+1; 
		        c<-0;
		        for (j in 1:length(s.c)){ 
			        if (is.na(s.c.begin[j]) || is.na(s.c.end[j])){
			        }else{
				        if (s.c.begin[j] <= ii && ii <= s.c.end[j]){
					        c<-c+1
				        }
			        } 
		        }
		        tt[t]<-ii
		        jj[t]<-c
	        }
	        lines(tt,jj,col=colormap[i],lwd=3)
            if (plotConcurrentMax == TRUE){
                points(xx<-tt[jj == max(jj)], yy<-jj[jj == max(jj)], cex=0.5, pch=22)
                text(xx<-tt[jj == max(jj)], yy<-jj[jj == max(jj)], max(jj), cex=0.5, pos=3)
            }
        }
    }
}
