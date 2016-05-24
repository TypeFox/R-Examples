cmp.rel <-  function (formula = formula(data), data = parent.frame(), ratetable = relsurv::slopop, 
     na.action,tau,conf.int=0.95,precision=1,add.times) 
    
    #formula: for example Surv(time,cens)~1 #not implemented for subgroups - DO IT!
    #data: the observed data set
    #ratetable: the population mortality tables
    #conf.type: confidence interval calculation (plain, log or log-log)
    #conf.int: confidence interval 
    #tau: max. cas do katerega racuna
{

    call <- match.call()
    rform <- rformulate(formula, data, ratetable, na.action)			#get the data ready
    data <- rform$data								#the data set
    se.fac <- sqrt(qchisq(conf.int, 1))						#factor needed for confidence interval

    if(missing(tau)) tau<-max(rform$Y)

    p <- rform$m								#number of covariates
    if (p > 0) #if covariates
    data$Xs <- strata(rform$X[, ,drop=FALSE ])					#make strata according to covariates
    else data$Xs <- rep(1, nrow(data))						#if no covariates, just put 1
    
    tab.strata <- table(data$Xs)						#unique strata values
    ntab.strata <- length(tab.strata)						#number of strata
    
    dtemp <- list(NULL)
    out <- as.list(rep(dtemp,ntab.strata*2))

    for (kt in 1:ntab.strata) {							#for each stratum
        inx <- which(data$Xs == names(tab.strata)[kt])				#individuals within this stratum
	
	extra <- as.numeric(seq(1,max(rform$Y[inx]),by=precision))
	if(!missing(add.times)) extra <- c(extra,as.numeric(add.times))
	tis <- sort(unique(pmin(tau,union(rform$Y[inx],extra))) )	#1-day long intervals used - to take into the account the continuity of the pop. part

	#if(!all.times)tis <- sort(unique(pmin(rform$Y[inx],tau)))					#unique times
	#else{
	#	tis <- sort(union(rform$Y[inx], as.numeric(1:floor(max(rform$Y[inx]))))) 	#1-day long intervals used - to take into the account the continuity of the pop. part
	#	tis <- unique(pmin(tis,tau))	
   	#}
	k <- length(tis)
	
	out[[2*kt-1]]$time <- out[[2*kt]]$time <- c(0,tis)
	
	temp <- exp.prep(rform$R[inx,,drop=FALSE],rform$Y[inx],ratetable,rform$status[inx],times=tis,fast=TRUE,cmp=T)		#calculate the values for each interval of time

	areae <- sum(temp$areae)/365.241		# sum(diff(c(0,tis))*temp$cumince)/365.241
	areap <- sum(temp$areap)/365.241		#sum(diff(c(0,tis))*temp$cumincp)/365.241
  
  	options(warn=-1)
  	out[[2*kt-1]]$est <- c(0,temp$cumince)
  	out[[2*kt-1]]$var <- c(0,temp$ve)
  	out[[2*kt-1]]$lower <- temp$cumince-se.fac*sqrt(temp$ve)
  	out[[2*kt-1]]$upper <- temp$cumince+se.fac*sqrt(temp$ve)
  	out[[2*kt-1]]$area <- areae
  	
  	out[[2*kt]]$est <- c(0,temp$cumincp)
	out[[2*kt]]$var <- c(0,temp$vp)
	out[[2*kt]]$lower <- temp$cumincp-se.fac*sqrt(temp$vp)
	out[[2*kt]]$upper <- temp$cumincp+se.fac*sqrt(temp$vp)
	out[[2*kt]]$area <- areap
	options(warn=0)
	
	ne <- sum(temp$ve<0)
	if(ne>0) warning(paste(names(tab.strata)[kt],": The estimated variance of crude mortality is negative in ", ne, " out of ", length(temp$ve)," intervals"), call. = FALSE)
	
	if(!missing(add.times)){
	out[[2*kt-1]]$index <- out[[2*kt]]$index <- unique(c(1,which(tis %in% c(rform$Y[inx],add.times,tau))))
	out[[2*kt-1]]$add.times <- out[[2*kt]]$add.times <- add.times
	}
	else out[[2*kt-1]]$index <- out[[2*kt]]$index <- unique(c(1,which(tis %in% c(rform$Y[inx],tau))))
	
   }
   if(p>0)names(out) <- paste(rep(c("causeSpec","population"),ntab.strata),rep(names(tab.strata),each=2))
   else names(out) <- c("causeSpec","population")
   class(out) <- "cmp.rel"
   out
}


plot.cmp.rel <- function (x, main = " ", curvlab, ylim = c(0, 1), xlim, wh = 2, 
    xlab = "Time (days)", ylab = "Probability", lty = 1:length(x), xscale=1,
    col = 1, lwd = par("lwd"), curves, conf.int, all.times=FALSE,...) 
{
#wh= upper left coordinates of the legend, if of length 1, the legend is placed top left.
    nc <- length(x)			#number of curves
    if (length(lty) < nc)		#if not enough different line types 
        lty <- rep(lty[1], nc)
    else lty <- lty[1:nc]
    if (length(lwd) < nc)		#if not enough different line widths 
        lwd <- rep(lwd[1], nc)
    else lwd <- lwd[1:nc]
    if (length(col) < nc)		#if not enough different colors 
        col <- rep(col[1], nc)
    else col <- col[1:nc]
    if (missing(curvlab)) {		#if no curve labels desired
        if (mode(names(x)) == "NULL") { #and no curve labels prespecified
            curvlab <- as.character(1:nc)
        }
        else curvlab <- names(x)[1:nc]	#use prespecified if they exist
    }
    if (missing(xlim)) {		#if no limits desired
        xmax <- 0
        for (i in 1:nc) {
            xmax <- max(c(xmax, x[[i]][[1]]/xscale))	#take max time over all strata
        }
        xlim <- c(0, xmax)
    }
    if(all.times){
    for(it in 1:nc){
       x[[it]]$index <- 1:length(x[[it]][[1]])
    }
    }
    if(missing(curves))curves <- 1:nc
    if(missing(conf.int))conf.int <- NULL
    curves <- unique(curves)
    conf.int <- unique(conf.int)
    if(any((curves %in% 1:nc)==FALSE)) stop(paste("The curves argument should be specified as a vector of integers from 1 to", nc,sep=" "))
    if(any((conf.int %in% 1:nc)==FALSE)) stop(paste("The conf.int argument should be specified as a vector of integers from 1 to", nc,sep=" "))
    if(any((conf.int %in% curves)==FALSE)) stop("Confidence interval may only be plotted if the curve is plotted, see argument curves")
    
    col_nums <- floor(seq(from=95,to=50,length.out=length(conf.int)+2))
    col.conf.temp <-  sapply(col_nums,function(x)paste("gray",as.character(x),sep=""))
    col.conf.int <- rep("white",nc)
    col.conf.int[conf.int] <-  col.conf.temp[-c(1,length(col.conf.temp))]

    plot((x[[1]][[1]]/xscale)[x[[1]]$index], (x[[1]][[2]])[x[[1]]$index], type = "n", ylim = ylim, xlim = xlim, 
        main = main, xlab = xlab, ylab = ylab, bty = "l", ...)		#plot estimates [[1]]=time, [[2]]=est
    if (length(wh) != 2) {
        wh <- c(xlim[1], ylim[2])
    }
    u <- list(...)
    if (length(u) > 0) {
        i <- pmatch(names(u), names(formals(legend)), 0)
        do.call("legend", c(list(x = wh[1], y = wh[2], legend = curvlab[curves], 
            col = col[curves], lty = lty[curves], lwd = lwd[curves], bty = "n", bg = -999999), 
            u[i > 0]))
    }
    else {
        do.call("legend", list(x = wh[1], y = wh[2], legend = curvlab[curves], 
            col = col[curves], lty = lty[curves], lwd = lwd[curves], bty = "n", bg = -999999))
    }
    for(i in conf.int){
    	 if(i%%2==0)with(x[[i]],polygon(c(time[index][!is.na(lower[index])],rev(time[index][!is.na(upper[index])]))/xscale,c(lower[index][!is.na(lower[index])],rev(upper[index][!is.na(upper[index])])),col = col.conf.int[i] , border = FALSE))
    	 else with(x[[i]],my.poly(time[index][!is.na(lower[index])]/xscale,time[index][!is.na(upper[index])]/xscale,lower[index][!is.na(lower[index])],upper[index][!is.na(upper[index])],col = col.conf.int[i] , border = FALSE))
    }
    for (i in curves) {
        tip <- "s"
        if(i%%2==0)tip <- "l"
        lines((x[[i]][[1]]/xscale)[x[[i]]$index], (x[[i]][[2]])[x[[i]]$index], lty = lty[i], col = col[i], 
            lwd = lwd[i], type=tip, ...)
    }
}

my.poly <- function(x1,x2,y1,y2,...){
	x1 <- rep(x1,each=2)[-1]
	y1 <- rep(y1,each=2)[-(2*length(y1))]
	x2 <- rep(x2,each=2)[-1]
	y2 <- rep(y2,each=2)[-(2*length(y2))]
	polygon(c(x1,rev(x2)),c(y1,rev(y2)),...)
}


print.cmp.rel <- function (x, ntp = 4, maxtime,xscale=365.241, ...) 
{
    nc <- length(x)

    if (missing(maxtime)) {
        maxtime <- 0
        for (i in 1:nc) maxtime <- max(maxtime, x[[i]]$time)
    }
    tp <- pretty(c(0, maxtime/xscale), ntp + 1)
    tp <- tp[-c(1, length(tp))]
    
    if(length(x[[1]]$add.times)>0 & length(x[[1]]$add.times)<5){
	tp <- sort(unique(c(tp,round(x[[1]]$add.times/xscale,1))))
    }
    cat("Estimates, variances and area under the curves:\n")
    print(checktimes(x, tp,xscale,area=TRUE), ...)
    invisible()
}

checktimes <- function (w, times,xscale=1,area=FALSE) 
{
    ng <- length(w)
    times <- sort(unique(times))*xscale
    nt <- length(times)
    storage.mode(times) <- "double"
    storage.mode(nt) <- "integer"
    ind <- matrix(0, ncol = nt, nrow = ng)
    oute <- matrix(NA, ncol = nt, nrow = ng)
    outv <- oute
    outa <- matrix(NA,ncol=1,nrow=ng)
    storage.mode(ind) <- "integer"
    slct <- rep(TRUE, ng)
    for (i in 1:ng) {
        if (is.null((w[[i]])$est)) {
            slct[i] <- FALSE
        }
        else {
            z <- rep(NA,nt)
            for(kt in 1:nt)z[kt] <- rev(which(w[[i]][[1]]<=times[kt]))[1]
            ind[i, ] <- z
            oute[i, ind[i, ] > 0] <- w[[i]][[2]][z]
            outa[i,] <- w[[i]][[6]]
            if (length(w[[i]]) > 2) 
                outv[i, ind[i, ] > 0] <- w[[i]][[3]][z]
        }
    }
    dimnames(oute) <- list(names(w)[1:ng], as.character(times/xscale))
    dimnames(outv) <- dimnames(oute) 
    rownames(outa) <- rownames(oute)
    colnames(outa) <- paste("Area at tau") 
    if(area)list(est = oute[slct, , drop = FALSE], var = outv[slct, , 
        drop = FALSE], area=outa[slct,,drop=FALSE])
    else list(est = oute[slct, , drop = FALSE], var = outv[slct, , 
        drop = FALSE])
}