"XvsYbarplot"=function (var1, var2, dataset, width = 1, space = NULL, names.arg = NULL, 
    legend.text = NULL, horiz = FALSE, density = NULL, 
    angle = 45, col = NULL, border = par("fg"), main = NULL, 
    sub = NULL, xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, 
    xpd = TRUE, log = "", axes = TRUE, axisnames = TRUE, cex.axis = par("cex.axis"), 
    cex.names = par("cex.axis"), inside = TRUE, plot = TRUE, 
    axis.lty = 0, offset = 0, add = FALSE, ...) 
{		
	dataset<- as.data.frame(dataset)
	if (is.character(var1)& is.character(var2)){
		num_var1 <- match(var1,names(dataset))
		num_var2 <- match(var2,names(dataset))
		height<-table(dataset[,num_var1],dataset[,num_var2])
		}
	else { height <- table(var1,var2)}

    if (!missing(inside)) 
        .NotYetUsed("inside", error = FALSE)
    if (is.null(space)) 
        space <- if (is.matrix(height)) 
            c(0, 1)
        else 0.2
    space <- space * mean(width)
    if (plot && axisnames && is.null(names.arg)) 
        names.arg <- if (is.matrix(height)) 
            colnames(height)
        else names(height)
	
    if (is.vector(height) || (is.array(height) && (length(dim(height)) == 
        1))) {
        height <- cbind(height)
               if (is.null(col)) 
            if (is.null(col))
		col<-rep(0)			
		NR<-length(height)		
		max<-max(height)		
		for (i in 1:NR){
			if (height[i]==max){ coli<-hsv(h=1,s=1,v=1,1,1)
						   col<-cbind(col,coli)}
			else{
					coli<-hsv(h=1,s=0.4+(height[i]/max)*0.6,v=1,1,1)	
					col<-c(col,coli)}
				}						 
		col<-col[-1]		
 }
    else if (is.matrix(height)) {
        if (is.null(col)) 		
		NR <- nrow(height)	
		NC <- ncol(height)            
		col=rep(0)
		a<-1/NC			
			
			for(i in 1:NC){	
						max<-max(height[,i])
						coli<-rep(0,NR)			
							for (j in 1:NR){
								if (height[j,i]==max){
												coli[j]<-hsv(h=a*i,s=1,v=1,1,1)
												
											    }
						else {	
										coli[j]<-hsv(h=a*i,s=0.4+(height[j,i]/max)*0.6,v=1,1)	
							}
									    }
				col<-c(col,coli)
					  }
			col<-col[-1]					
    }

    else stop("'height' must be a vector or a matrix")
    if (is.logical(legend.text)) 
        legend.text <- if (legend.text && is.matrix(height)) 
            rownames(height)
    stopifnot(is.character(log))
    logx <- logy <- FALSE
    if (log != "") {
        logx <- length(grep("x", log)) > 0L
        logy <- length(grep("y", log)) > 0L
    }
    if ((logx || logy) && !is.null(density)) 
        stop("Cannot use shading lines in bars when log scale is used")
    NR <- nrow(height)
    NC <- ncol(height)
            if (length(space) == 2) 
            space <- rep.int(c(space[2], rep.int(space[1], NR - 
                1)), NC)
        width <- rep(width, length.out = NR)
        offset <- rep(as.vector(offset), length.out = length(width))
    delta <- width/2
    w.r <- cumsum(space + width)
    w.m <- w.r - delta
    w.l <- w.m - delta
    log.dat <- (logx && horiz) || (logy && !horiz)
    if (log.dat) {
        if (min(height + offset, na.rm = TRUE) <= 0) 
            stop("log scale error: at least one 'height + offset' value <= 0")
        if (logx && !is.null(xlim) && min(xlim) <= 0) 
            stop("log scale error: 'xlim' <= 0")
        if (logy && !is.null(ylim) && min(ylim) <= 0) 
            stop("log scale error: 'ylim' <= 0")
        rectbase <- if (logy && !horiz && !is.null(ylim)) 
            ylim[1]
        else if (logx && horiz && !is.null(xlim)) 
            xlim[1]
        else 0.9 * min(height, na.rm = TRUE)
    }
    else rectbase <- 0
    
    rAdj <- offset + (if (log.dat) 
        0.9 * height
    else -0.01 * height)
    delta <- width/2
    w.r <- cumsum(space + width)
    w.m <- w.r - delta
    w.l <- w.m - delta
    
    num_mod <- nlevels(var1)
    if (horiz) {
        if (is.null(xlim)) 
            xlim <- range(rAdj, height + offset, na.rm = TRUE)
        if (is.null(ylim)) 
            ylim <- c(min(w.l), max(w.r) + num_mod + (num_mod-1))
    }
    else {
        if (is.null(xlim)) 
            xlim <- c(min(w.l), max(w.r))
        if (is.null(ylim)) 
            ylim <- range(rAdj, height + offset + num_mod*5, na.rm = TRUE)
    }
            w.m <- matrix(w.m, ncol = NC)
	par(mar = par("mar") + c(1,0,0,0))
    if (plot) {
        opar <- if (horiz) 
            par(xaxs = "i", xpd = xpd)
        else par(yaxs = "i", xpd = xpd)
        on.exit(par(opar))
        if (!add) {
            plot.new()
		if(horiz){
			if (is.character(attributes(var1)$levels)==TRUE){
			if (max(nchar(attributes(var1)$levels))>8){
			par(mar = par("mar") + c(0,round(max(nchar(attributes(var1)$levels))/3),0,0))}
											}
			else {par(mar=par("mar")+c(0,5,0,0))}					
			plot.window(xlim, ylim, log = log, ...)
				}
            else  plot.window(xlim, ylim, log = log, ...)
        }
        xyrect <- function(x1, y1, x2, y2, horizontal = TRUE, 
            ...) {
            if (horizontal) 
                rect(x1, y1, x2, y2, ...)
            else rect(y1, x1, y2, x2, ...)
        }
        	
            xyrect(rectbase + offset, w.l, c(height) + offset, 
                w.r, horizontal = horiz, angle = angle, density = density, 
                col = col, border = border)
	      
        if (axisnames && !is.null(names.arg)) {
            at.l <- if (length(names.arg) != length(w.m)) {
                if (length(names.arg) == NC) 
                  colMeans(w.m)
                else stop("incorrect number of names")
            }
            else w.m
		if (!horiz)
	if (is.character(attributes(var1)$levels)==TRUE){
		for (i in 1:nlevels(var2)){
			if (nchar(attributes(var2)$levels[i]) > 11) names.arg[i] <- substring(attributes(var2)$levels[i],1,11)
						}								
			}
		if(horiz){
			axis(2, at = at.l, labels = names.arg, lty=axis.lty,cex.axis = cex.names, las = 2)}
		else {
			axis(1, at = at.l, labels = names.arg, lty=axis.lty,cex.axis = cex.names, las = 0)}
                  }
        if (!is.null(legend.text)) {
            legend.col <- rep(col, length.out = length(legend.text))
            if (!horiz)  {
		    legend.text <- legend.text
                legend.col <- legend.col
                density <- rev(density)
                angle <- rev(angle)  		
		}
		num.legend <- c("1st","2nd","3rd")
		for (i in 4:20){
			num.legendi <- c(paste(i,"th"))
			num.legend <- c(num.legend,num.legendi)}
		num.legend <- num.legend[1:dim(height)[1]]
            xy <- par("usr")
            if (horiz){
		legend2( xy[2] - xinch(0.1), xy[4] - yinch(0.1), legend = paste(num.legend,"bar:",legend.text), 
                angle = angle, density = density, fill = legend.col, bty="n", cex = 1 - 0.04*num_mod, 
                xjust = 1, yjust = 1)}
		else {
		legend2( xy[2] - xinch(0.1), xy[4] - yinch(0.1), legend = paste(num.legend,"bar:",legend.text), 
                angle = angle, density = density, fill = legend.col, bty="n",
                xjust = 1, yjust = 1)}
		        }
        if (is.character(var1)& is.character(var2)){
        title(main = paste(names(dataset[num_var1]),"depending on", names(dataset[num_var2])), sub = sub, xlab = xlab, ylab = ylab, 
            ...)}
	else {
		for (i in 1:length(dataset)){
			a <- match(var1,dataset[,i])
			b <- match(var2,dataset[,i])
				if (any(is.na(a)) == FALSE){ rep1=i}	
				if (any(is.na(b)) == FALSE){ rep2=i}		
					 		 }
					title(main = paste(names(dataset[rep1]),"depending on", names(dataset[rep2])), sub = sub, xlab = xlab, ylab = ylab, 
            ...)}
	
        if (axes)
		if (horiz){
			axis(1, cex.axis = cex.axis, las = 0, ...)}
		else { axis(2, cex.axis = cex.axis, las = 2, ...)}
                   invisible(w.m)
    }
    else w.m

}





#barplot2(var1="Sport",var2="age_Q",tea,legend.text=TRUE,horiz=TRUE)
#X11()
#barplot2(var1=tea[,18],var2=tea[,20],tea,legend.text=TRUE,cex.names=1,horiz=TRUE)
