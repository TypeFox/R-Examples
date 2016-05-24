

#-----------------
#A function to spread out the data. Also returns the min and max of the spread out data, for each rater and level.


factor_spreadout <- function(d){
	d <- as.data.frame(d)
	n <- nrow(d)
	dint <-apply(d,2,rank,ties.method="first") #maps the categorical variables to 1..n
	dint <- dint/n
	dtab <- lapply(d, table)
	tops <- lapply(dtab, function(x) cumsum(x)/n)
	bots <- lapply(tops, function(x) c(0,x[-length(x)]))
	bars <- mapply(cbind,bots,tops,SIMPLIFY=FALSE)
    for (i in 1:length(tops))
      rownames(bars[[i]]) <- names(tops[[i]])
	return(list(data=dint,bars=bars))
	}
	


rater_spreadout <- function(d,levs, minspace=NULL, scale=FALSE){
	d <- as.data.frame(d)
	dtab <- sapply(d, function(x) table(factor(x,levels=levs)))

   if (is.null(minspace))
	minspace <- max((dtab[-1,]+dtab[-nrow(dtab),])/2)*1.1
	#minspace  <- max(dtab)
     barc <- row(dtab)    
     barw <- dtab/minspace
    barb <- barc - barw/2
    bart <- barc+barw/2
    if (scale){
    	m1 <- rep(min(barb),ncol(barb))
    	m2 <- rep(max(bart),ncol(bart))
    	barb <- scale(barb,center=m1,scale=m2)
    	bart <- scale(bart,center=m1,scale=m2)
    	barc <- scale(barc,center=m1,scale=m2)
    	}
	dnew <- d
	for (j in 1:ncol(dtab)){
		  x <- d[,j]
	      for (i in 1:nrow(dtab)){
			n <- dtab[i,j]
			k <- levs[i]
			spread <- seq(barb[i,j], bart[i,j],length.out= n)
			if (n > 1)
		    dnew[x==k,j] <- spread		
		    else if (n==1)  dnew[x==k,j] <- barc[i,k] }
		}
	bars <- list(NULL)
    for (i in 1:ncol(barb))
       bars[[i]] <- cbind(barb[,i], bart[,i])
	return(list(data=dnew,bars=bars))
	}
	


axis_bars <- function(bars,o,barvars=sort(unique(o)),width=.2){
	barvars.o <- o[o %in% barvars]
	barvars.x <- seq(along=o)[o %in% barvars]
	nbars <- sapply(bars, nrow)
    left <- rep(barvars.x,times=nbars[barvars.o])
    
    left <- left - width/2
    right <- left + width
    bars.o <- bars[barvars.o]
    barnames <- unlist(lapply(bars.o,rownames))
    barb <- unlist(lapply(bars.o, function(b) b[,1]))
    bart <- unlist(lapply(bars.o, function(b) b[,2]))
	ans <- cbind(left,barb, right, bart)
	rownames(ans) <- barnames
	colnames(ans) <- c("left","bottom", "right","top")
	ans
	}
	
	


catpcp <- function (data, order = NULL, pcpbars, barvars=1:ncol(data), pcpbars.border="black",pcpbars.col=NULL,pcpbars.labels=FALSE,pcpbars.axis.at=NULL,pcpbars.axis.labels=NULL,axis.width=.2,connect=TRUE,...) {
    	
    	pcp(data,order,axis.width=axis.width,connect=connect,...)
    	oldxpd <- par("xpd")
        par("xpd"=TRUE) 
        if (is.null(order)) 
          order <- 1:ncol(data)
        else
          data <- data[, order]
        
        pcpbars <- axis_bars(pcpbars,order,barvars,axis.width)


        rect(pcpbars[,1],pcpbars[,2],pcpbars[,3],pcpbars[,4],col=pcpbars.col,
            border=pcpbars.border)
        if (pcpbars.labels)
        text((pcpbars[,1]+pcpbars[,3])/2,(pcpbars[,2]+pcpbars[,4])/2,rownames(pcpbars))
       if (!is.null(pcpbars.axis.at))
           axis(2,las=2,lwd=0,labels=pcpbars.axis.labels,at=pcpbars.axis.at)
        par("xpd"=oldxpd)
     }
