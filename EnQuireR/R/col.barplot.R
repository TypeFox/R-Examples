"col.barplot"=function (height, width = 1, space = NULL, names.arg = NULL,
    legend.text = NULL, beside=TRUE, horiz = FALSE, density = NULL, 
    angle = 45, col = NULL, border = par("fg"), main = NULL, 
    sub = NULL, xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, 
    xpd = TRUE, log = "", axes = TRUE, axisnames = TRUE, cex.axis = par("cex.axis"), 
    cex.names = par("cex.axis"), inside = TRUE, plot = TRUE, 
    axis.lty = 0, offset = 0, add = FALSE, ...) 
{
    if (!missing(inside)) 
        .NotYetUsed("inside", error = FALSE)
    if (is.null(space)) 
        space <- if (is.matrix(height) ) 
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
        beside <- TRUE 			
	  if (is.null(col))
		col<-rep(0)			## on crée un répertoire de 0
		NR<-length(height)	## on veut le nombre de lignes maximum	
		max<-max(height)		## on veut le maximum parmi les modalités de la variable
		for (i in 1:NR){
			if (height[i]==max){ coli<-hsv(h=1,s=1,v=1,1,1)
						   col<-cbind(col,coli)}
			else{
					coli<-hsv(h=1,s=0.4+(height[i]/max)*0.6,v=1,1,1)	## pour chaque modalité de la variable, on affecte une teinte de couleur proportionnelle à sa représentativité par rapport à la modalité qui a l'effectif maximumu.
					col<-c(col,coli)}
				}						## on met dans le vecteur col les couleurs pour chaque modalité. 
		col<-col[-1]		##on retire la première coordonnée qui contient un zéro.
 }
    else if (is.matrix(height)) {
        if (is.null(col)) 		
		NR <- nrow(height)	## on veut le nombre de lignes	
		NC <- ncol(height)      ## on veut le nombre de colonnes      
		col=rep(0)
		a<-1/NC			## permet de diviser la palette de couleurs en autant de couleurs qu'on a de colonnes.
			
			for(i in 1:NC){	
						max<-max(height[,i])
						coli<-rep(0,NR)			## pour chaque colonne de la matrice, on crée un répertoire qui va répertorier les différentes teintes dégradées de la couleur de la colonne i.
							for (j in 1:NR){
								if (height[j,i]==max){
												coli[j]<-hsv(h=a*i,s=1,v=1,1,1)
												
											    }
						else {	
										coli[j]<-hsv(h=a*i,s=0.4+(height[j,i]/max)*0.6,v=1,1)	## pour toutes les colonnes, on fait comme expliqué 13 lignes plus haut
							}
									    }
				col<-c(col,coli)
					  }
			col<-col[-1]		## on retire la première coordonnée qui vaut 0.			
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
    if (beside) {
        if (length(space) == 2) 
            space <- rep.int(c(space[2], rep.int(space[1], NR - 
                1)), NC)		
        width <- rep(width, length.out = NR)	
    }
    else {
        width <- rep(width, length.out = NC)
    }
    offset <- rep(as.vector(offset), length.out = length(width))	
    delta <- width/2
    w.r <- cumsum(space + width)	##donne la longueur totale
    w.m <- w.r - delta
    w.l <- w.m - delta			##donne le nombre de variables.
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
    if (!beside) 
        height <- rbind(rectbase, apply(height, 2, cumsum))
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
            ylim <- c(min(w.l), max(w.r))
    }
    else {
        if (is.null(xlim)) 
            xlim <- c(min(w.l), max(w.r))
        if (is.null(ylim)) 
            ylim <- range(rAdj, height + offset, na.rm = TRUE)
    }
    if (beside) 
        w.m <- matrix(w.m, ncol = NC)
    if (plot) {
        opar <- if (horiz) 
            par(xaxs = "i", xpd = xpd)
        else par(yaxs = "i", xpd = xpd)
        on.exit(par(opar))
        if (!add) {
            plot.new()
            plot.window(xlim, ylim, log = log, ...)
        }
        xyrect <- function(x1, y1, x2, y2, horizontal = TRUE, 
            ...) {
            if (horizontal) 
                rect(x1, y1, x2, y2, ...)
            else rect(y1, x1, y2, x2, ...)
        }
        if (beside) 
            xyrect(rectbase + offset, w.l, c(height)+ offset, 
                w.r, horizontal = horiz, angle = angle, density = density, 
                col = col, border = border)					##trace les bares quand beside=TRUE
   			
	  else {
            for (i in 1:NC) {
                xyrect(height[1:NR, i] + offset[i], w.l[i], height[-1, 
                  i] + offset[i], w.r[i], horizontal = horiz, 
                  angle = angle, density = density, col = col, 
                  border = border)							##trace les bares quand beside=FALSE
            		    }
        	}
        if (axisnames && !is.null(names.arg)) {
            at.l <- if (length(names.arg) != length(w.m)) {
                if (length(names.arg) == NC) 
                  colMeans(w.m)
                else stop("incorrect number of names")
            }
            else w.m
            axis(if (horiz) 
                2
            else 1, at = at.l, labels = names.arg, lty = axis.lty, 
                cex.axis = cex.names, ...)
        }
        if (!is.null(legend.text)) {
            legend.col <- rep(col, length.out = length(legend.text))
            if ((horiz ) || (!horiz))  {
		    legend.text <- rev(legend.text)
                legend.col <- rev(legend.col)
                density <- rev(density)
                angle <- rev(angle)
            }
            xy <- par("usr")
            legend(xy[2] - xinch(0.1), xy[4] - yinch(0.1), legend = legend.text, 
                angle = angle, density = density, fill = legend.col, 
                xjust = 1, yjust = 1)
        }
        title(main = main, sub = sub, xlab = xlab, ylab = ylab, 
            ...)
        if (axes) 
            axis(if (horiz) 
                1
            else 2, cex.axis = cex.axis, ...)
        invisible(w.m)
    }
    else w.m
}
