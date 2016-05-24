"barplotclustervsvar"=function (var, cluster, dataset, width = 1, space = NULL, names.arg = NULL, 
     horiz = FALSE, density = NULL, 
    angle = 45, border = par("fg"), 
    sub = NULL, xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, 
    xpd = TRUE, log = "", axes = TRUE, cex.axis = par("cex.axis"), 
    cex.names = par("cex.axis"), inside = TRUE, plot = TRUE, 
    axis.lty = 0, offset = 0, add = FALSE, ...) 
{
	dataset <- as.data.frame(dataset)
	if (is.character(var)){
		num_var <- match(var,names(dataset))
		height <- table(dataset[,num_var],cluster)
		}
	else {height <- table(var,cluster)}
	
    if (!missing(inside)) 
        .NotYetUsed("inside", error = FALSE)
    if (is.null(space)) 
        space <- if (is.matrix(height) ) #&& beside) 
            c(0, 1)
        else 0.2
    space <- space * mean(width)
    if (plot && is.null(names.arg)) 
        names.arg <- if (is.matrix(height)) 
            colnames(height)
        else names(height)

        
            col <- col <- c("black", "red", "green3", "blue", "cyan", "magenta",
            "darkgray", "darkgoldenrod", "darkgreen", "violet",
            "turquoise", "orange", "lightpink", "lavender", "yellow",
            "lightgreen", "lightgrey", "lightblue", "darkkhaki",
            "darkmagenta", "darkolivegreen", "lightcyan", "darkorange",
            "darkorchid", "darkred", "darksalmon", "darkseagreen",
            "darkslateblue", "darkslategray", "darkslategrey",
            "darkturquoise", "darkviolet", "lightgray", "lightsalmon",
            "lightyellow", "maroon")
		
		col2 <- matrix(0,3,length(col))
		col2 <- col2rgb(col)
NR <- nrow(height)
NC <- ncol(height)
grey.colors <- grey(seq(0,0.7,0.05))
col_print <- rep(0)
max1 <- max(height[,1])
min1 <- min(height[,1]) 
h <- (max1-min1)/15
col_print1 <- rep(0,nlevels(var))
  for (k  in 1:NR){	
		if (height[k,1]==max1){ 
						col_print1[k] <- grey.colors[1]}
		else { az <- max1-height[k,1]
			 er <- floor(az/h)
				if(er == 0) {col_print1[k]<- grey.colors[2]}			 
			      else {col_print1[k] <- grey.colors[er]}
			} 
					}
col_print <- c(col_print,col_print1)


	for (i in 2:NC){
		max <- max(height[,i])
		col3i <- matrix(0,3,nlevels(var))
				col_printi <- rep(0,nlevels(var))		
			for (j in 1:NR){
					if (height[j,i]==max){ col3i[,j] <- rgb2hsv(col2[,i])
									col_printi[j] <-  hsv( h = col3i[1,j], s = col3i[2,j], v = col3i[3,j], 1)
									}
					else { col3i[,j] <- rgb2hsv(col2[,i])
						 col_printi[j] <- hsv(h = col3i[1,j], s = 0.4+(height[j,i]/max)*0.6,v = 1, 1) 
								}
						}
				col_print <- c(col_print,col_printi)}
		col_print <-col_print[-1]
	
        
        legend.text <- rownames(height)					 
            
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

    if (horiz) {
        if (is.null(xlim)) 
            xlim <- range(rAdj, height + offset, na.rm = TRUE)
        if (is.null(ylim)) 
            ylim <- c(min(w.l), max(w.r)+ NR + (NR-1))
    }
    else {
        if (is.null(xlim)) 
            xlim <- c(min(w.l), max(w.r))
        if (is.null(ylim)) 
            ylim <- range(rAdj, height + offset + NR*7, na.rm = TRUE)
    }
     
        w.m <- matrix(w.m, ncol = NC)
    if (plot) {
        opar <- if (horiz) 
            par(xaxs = "i", xpd = xpd)
        else par(yaxs = "i", xpd = xpd)
        on.exit(par(opar))
        if (!add) {
            plot.new()
            plot.window(xlim, ylim, log = log)
        }
        xyrect <- function(x1, y1, x2, y2, horizontal = TRUE, 
            ...) {
            if (horizontal) 
                rect(x1, y1, x2, y2, ...)
            else rect(y1, x1, y2, x2, ...)
        }
         
            xyrect(rectbase + offset, w.l, c(height) + offset, 
                w.r, horizontal = horiz, angle = angle, density = density, 
                col = col_print, border = border)
               }
        if (!is.null(names.arg)) {
            at.l <- if (length(names.arg) != length(w.m)) {
                if (length(names.arg) == NC) 
                  colMeans(w.m)
                else stop("incorrect number of names")
            }
            else w.m
            if (horiz) {   
            	axis(2, at = at.l, labels = paste("group",names.arg), lty = axis.lty, cex.axis = cex.names, las = 2)
				}
		else{
            	axis(1, at = at.l, labels = paste("group",names.arg), lty = axis.lty, cex.axis = cex.names, las = 0)
			} 
        }
        #if (!is.null(legend.text)) {
            legend.col <- rep(col, length.out = length(legend.text))
            if (!horiz) {
                legend.text <- legend.text
                legend.col <- legend.col
                density <- rev(density)
                angle <- rev(angle)
            }
		else { legend.text <- rev(legend.text)}

      num.legend <- c("1st","2nd","3rd")
		for (i in 4:20){
			num.legendi <- c(paste(i,"th"))
			num.legend <- c(num.legend,num.legendi)}
		num.legend <- num.legend[1:dim(height)[1]]
	xy <- par("usr")
	if (horiz){
		legend2( xy[2] - xinch(0.1), xy[4] - yinch(0.1), legend = paste(num.legend,"bar:",legend.text), 
                angle = angle, density = density, fill = legend.col, bty="n", cex = 1 - 0.04*NR, 
                xjust = 1, yjust = 1)}
		else {
		legend2( xy[2] - xinch(0.1), xy[4] - yinch(0.1), legend = paste(num.legend,"bar:",legend.text), 
                angle = angle, density = density, fill = legend.col, bty="n",
                xjust = 1, yjust = 1)}
		       # }
	if (is.character(var)){ 
        title(main = paste(names(dataset[var]),"by cluster"), sub = sub, xlab = xlab, ylab = ylab, 
            ...)}
	else {	
			for (i in 1:length(dataset)){
					a <- match(var,dataset[,i]) 
				if (any(is.na(a)) == FALSE){ rep=i }
								}
		
	  title(main = paste(names(dataset[rep]),"by cluster"), sub = sub, xlab = xlab, ylab = ylab)} 
            

        if (axes) 
            if (horiz){axis(1, cex.axis = cex.axis, las = 0)}			
		else { axis(2, cex.axis = cex.axis, las = 2)}
        invisible(w.m)
    }
#    else w.m
#}

#barplot_cluster_vs_var(var=tea[,20],cluster=v,dataset=tea,legend.text=TRUE,horiz=F)
