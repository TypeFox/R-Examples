track <-
function(ephem, mag = 7, edge = 0.2, cex.star = 1,
         xlab = "Right Ascension", ylab = expression(paste("Declination", degree)), 
         col.track = "red", interval.lab = 15, lwd.track = 2, grid = TRUE, 
		 bright.lab = TRUE, bright, starcat, ...){
	
	
         background <- function(xlab = "",ylab = "", xlim = c(0, 24), 
         	                   ylim = c(-45, 45), xinterval = 2, yinterval = 15 , 
         	                   xunit = "h", yunit = expression(degree), axis = TRUE, 
         	                   auto = TRUE, grid = FALSE, ...){            
             par(xaxs = "i", yaxs = "i")
             ### The xticks of chart are always reversed.
             plot(0, 0, xlim = rev(xlim), ylim = ylim, axes = FALSE, xlab = xlab, ylab = ylab, type = "n")
         	
         		## Add axes for the current plot
             	xticks <- axTicks(1)
                  yticks <- axTicks(2)
             	
             	xlabs <- rev(paste(round((seq(0,24, by = abs(xticks[2]
                   		 - xticks[1]))),digits = 2), xunit, sep = ""))
                  xat   <- rev(round(seq(0,24, by = abs(xticks[2] - xticks[1])),2))
             	
             	ylabs <- paste(round(seq(-60,60, by = abs(yticks[2] 
         		         - yticks[1])), digits = 2), "", sep = "")
                  yat   <- round(seq(-60,60, by = abs(yticks[2] - yticks[1])),2)
             	
                  axis(1, at = xat, labels = xlabs)
                  axis(2, at = yat, labels = ylabs)
             
             if(grid){
                  abline(h = yat, col = "grey",  ...)
                  abline(v = xat, col = "grey",  ...)
             }
             ## the central line
             abline(v = 12)
             abline(h = 0)
             box()
         }
         		
         ### Add the tracks of a planet according to the input positions.
         tracks <- 
         function(position, limit = 1, col = 2, lwd = 2){
                     
             ###Find the intercept
             whichint <- 
             function(x, limit = 1){
                ind1 <- 1:(length(x)-1)
                ind2 <- ind1+1
                res <- x[ind1] - x[ind2]
                which(abs(res) > limit)
             }
             
             breaks <- whichint(x = position[,1], limit = 1)
             if(length(breaks) == 0){
             lines(position[,1],position[,2], col = col, lwd = lwd) 
         	## Lines have no breaks should not be cutted
             }else{
               for(i in 1:(length(breaks))){
                 if(i == 1){
                     brka <- 1
                     brkb <- breaks[i]   ## Beginning to the first cut.
                 }
                 if(i > 1){
                     brka <- breaks[i - 1] + 1
                     brkb <- breaks[i]   ## 
                 }
               lines(position[,1][brka:brkb],position[,2][brka:brkb], col = col, lwd = lwd)
              }
               brkc <- brkb +1
               brkd <- length(position[,1])
               lines(position[,1][brkc:brkd],position[,2][brkc:brkd], col = col, lwd = lwd)
             }
         }
             
         ## The magnitude of star
         magselect <- function(x, magnitude){
             res <- x[x[,3] < magnitude, ]
             res <- na.omit(res)
             return(res)
         }
             
         ### Add stars to the sky chart
         add.stars <- 
         function(catalog, xlim = c(0, 24), ylim = c(-30, 30), mag = 3, times = 1, ...){
             	## Subseting stars to be plot
             subres <- magselect(catalog, magnitude = mag)
             x1 <- range(xlim)
             x2 <- range(ylim)
             target <- subres[(subres[,1] >= x1[1])&(subres[,1] <= x1[2])&
         	                 (subres[,2] >= x2[1])&(subres[,2] <= x2[2]),]
             
             ## Calculating Size for each star
             if(is.integer(mag)){
                 magint <- 0:mag
             }else{
             magint <- c(0:floor(mag), mag)	 
             }
             int <- findInterval(target[,3], magint)
             ddd <- rep(0, length(int))
             pointsize <- (length(magint):1)/length(magint)
             for(i in 1:length(magint)){
               ind <- (magint[i] == int)
               ddd[ind] <- pointsize[i]
             }
             ## Add stars to a certain background
             points(x = (target[,1]), y = (target[,2]), xlim = xlim, ylim = ylim, 
         		       type = "p", cex = times * ddd, ...)
         }
             
             
         ### Add label for star tracks
         addlab.eqc <- 
         function(ephem, interval = 10, upper = 2, ...)
         {
             index <- seq(1, nrow(ephem), by = interval)
             x <- ephem[index,]
             points(x[,1], x[,2], pch = 19, col = 2, cex = 1.2, ...)
             text(x[,1], x[,2] + upper, substring(row.names(x), 6,10), ...)
         }
             
             
         ### Add Sun, Moon, and major planets to the chart. 
         add.planet <- function(time = jd(), col = 1, pcex = 3, pch = 21,...){
             labs <- c("Sun","Mon","Mercury","Venus","Mars","Jupiter",
         	          "Saturn","Uranus","Neptune","Pluto")
             position <- planets(time)
             text(position[,1], position[,2], labs, cex = 0.8)
             points(position[,1], position[,2],col = col, cex = pcex, pch = pch, ... )
         }
             
         #### Add ecliptic 
         add.ecliptic <- 
         function(){
             a <- ecc(1:360,rep(0,360),1:360)
             aa = as.eqc(a)
             lines((aa[,1])[-nrow(aa)], aa[,2][-(nrow(aa))])
         }
         
         ## Add labels of the bright stars
         add.brightlab <- 
         function(bright,upper){    
             starlab <- rownames(bright)
             text(bright[,1], bright[,2] + upper, starlab, cex = 0.7)
         }
         
         ### extend the width of xlim or ylim for plot
         egde.extent <- 
          function(x, edge  ) {
             x1 = x[1] - edge*(x[2]-x[1])
             x2 = x[2] + edge*(x[2]-x[1])
             x0 = c(x1, x2)
             return(x0)
         }
    
	xlim <- egde.extent(range(ephem[,1]), edge = edge)
	ylim <- egde.extent(range(ephem[,2]), edge = edge)
	
    if(xlim[1] < 0) xlim[1] =  0
    if(xlim[2] > 24) xlim[2] = 24
    
    if(xlim[1] < -90) xlim[1] = -90
    if(xlim[2] > 90) xlim[2] =   90

    par(tcl = 0.2)
	
	## plot the background
    background(xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, axis = axis, 
	          grid = grid,...)
	## add stars 		  
    add.stars(starcat, mag = mag, pch = 19, time = cex.star)
	
	if(bright.lab){
	add.brightlab(bright, upper = abs(ylim[2]-ylim[1])/20)
	}
	
    ## add tracks of planets
	tracks(ephem, col = col.track, lwd = lwd.track )
	## add labels of planets
    addlab.eqc(ephem, interval = interval.lab, upper = abs(ylim[2]-ylim[1])/20)
	
	nam <- rownames(ephem)
	astro <- substring(nam[1], 12, nchar(rownames(ephem[1,])))
	begin <- min(as.Date(substring(nam, 1,10))) 
	end <- max(as.Date(substring(nam, 1,10)))
	title(paste("Track of",astro, "from",begin,"to",end))
}

		        
