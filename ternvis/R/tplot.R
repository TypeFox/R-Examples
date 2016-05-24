tplot <-
function(x=cbind(1,1,1)/3, 
                    L = diag(c(1,1,1))/sqrt(2),    # score matrix is Brier by default
		    scale = 1, 
		    dimnames = NULL, 
		    dimnames_position = c("corner","none"), 
		    dimnames_color = "black", 
		    id = NULL, 
		    id_color = "black", 
                    coordinates = FALSE, 
		    grid = TRUE, 
		    grid_color = "gray", 
		    labels = c("inside","outside", "none"), 
		    labels_color = "darkgray", 
		    border = "grey", 
                    bg = "white", 
		    pch = 19, 
		    cex = 1, 
		    prop_size = FALSE, 
		    col = "red", 
                    main = "ternary plot", 
		    sub = NULL,
		    newpage = TRUE, 
		    pop = TRUE, 
		    col.main="black",
		    col.sub="black",
		    ...)
{
    oB <- rbind(1,0,0)
    oN <- rbind(0,1,0)
    oA <- rbind(0,0,1)
    a <- as.double(sqrt(t(oB-oN)%*%t(L)%*%L%*%(oB-oN)))
    b <- as.double(sqrt(t(oA-oN)%*%t(L)%*%L%*%(oA-oN)))
    n <- as.double(sqrt(t(oB-oA)%*%t(L)%*%L%*%(oB-oA)))
    phi <- acos((a^2+n^2-b^2)/(2*a*n))
    M32 <- rbind(cbind(0 ,a*cos(phi),n),
                 cbind(0 ,a*sin(phi),0)   )     
    colnames(M32) <- NULL

    labels <- match.arg(labels)
    
    if (grid == TRUE)
        grid <- "dotted"
	
    if (coordinates)
        id <- paste("(", 
	            round(x[,1] * scale, 1), 
		    ",", 
		    round(x[,2] * scale, 1),
		    ",", 
		    round(x[,3] * scale, 1),
		    ")",
                    sep = "")
	    
    dimnames_position <- match.arg(dimnames_position)
    
    if (is.null(dimnames) && dimnames_position != "none")
        dimnames <- colnames(x)
	
    if (is.logical(prop_size) && prop_size)
        prop_size <- 3
	
    if (ncol(x) != 3)
        stop("Need a matrix with 3 columns")
	
    if (any(x < 0))
       {
        print("X must be non-negative: resetting negative elements")
	x[x<0] <- 0
       }
	
    s <- rowSums(x)
    
    if (any(s <= 0))
        stop("STOP: each row of X must have a positive sum")
	
    x <- x/s
        
    if (newpage)
        grid.newpage()
	
    xmin <- min(0,a*cos(phi))
    xmax <- max(n,a*cos(phi))
    ymin <- 0
    ymax <- a*sin(phi)	
        
    height <- ymax - ymin
    width  <- xmax - xmin
    maxhw  <- max(height,width) 
    extra  <- 0.03*maxhw
    
    xlim <- c(xmin-extra, xmin + maxhw+extra)
    ylim <- c(ymin-extra, ymin + maxhw+extra)
    
    xoff <- 0
    yoff <- 0.5*(maxhw-height)
    
    pushViewport(viewport(width = unit(1, "snpc")))
    
    if (!is.null(main))
        grid.text(main, 
	          y = 0.9, 
		  gp = gpar(fontsize = 18, fontstyle = 1, col=col.main))

    if (!is.null(sub))
        grid.text(sub, 
	          y = 0.1, 
		  gp = gpar(fontsize = 18, fontstyle = 1, col=col.sub))
	
    pushViewport(viewport(width  = 0.8, 
                          height = 0.8, 
			  xscale = xlim,
                          yscale = ylim))
	
    eps <- 0.01
    top <- ymax
    
    grid.polygon(x = xoff + M32[1,], 
                 y = yoff + M32[2,], 
		 default.units="native",
		 gp = gpar(fill = bg, col = border), ...)

    if (dimnames_position == "corner") {
        grid.text(x = xoff + c(    0,  a*cos(phi),     n), 
	          y = yoff + c(-0.02,  top + 0.02, -0.02), 
		  label = dimnames, 
		  default.units="native",
		  gp = gpar(fontsize = 12))
    }






    if (is.character(grid))
        for (i in 1:4 * 0.2) {
            grid.lines(x=xoff + c(1 - i, (1 - i)/2),
	               y=yoff + c(0, 1 - i) * height,
                       gp = gpar(lty = grid, col = grid_color),
		       default.units="native")
            grid.lines(x=xoff + c(1 - i, 1 - i + i/2),
	               y=yoff + c(0, i) * height,
                       gp = gpar(lty = grid, col = grid_color),
		       default.units="native")
            grid.lines(x=xoff + c(i/2, 1 - i + i/2),
	               y=yoff + c(i, i) * height, 
		       gp = gpar(lty = grid,col = grid_color),
		       default.units="native")
		       
            if (labels == "inside") {
                grid.text(x = xoff + (1 - i) * 3/4 - eps, 
		          y = yoff + (1 - i)/2 * height, 
			  label = i * scale, 
			  gp = gpar(col = labels_color), 
			  rot = 120,
		       default.units="native")
                grid.text(x = xoff + 1 - i + i/4 + eps, 
		          y = yoff + i/2 * height - eps,
			  label = (1 - i) * scale, 
			  gp = gpar(col = labels_color),
                          rot = -120,
		       default.units="native")
                grid.text(x = xoff + 0.5, 
		          y = yoff + i * height + eps, 
			  label = i * scale, 
			  gp = gpar(col = labels_color),
		       default.units="native")
            }
            if (labels == "outside") {
                grid.text(x = xoff + (1 - i)/2 - 6 * eps, 
		          y = yoff + (1 - i) * height, 
			  label = (1 - i) * scale, 
			  gp = gpar(col = labels_color),
		       default.units="native")
                grid.text(x = xoff + 1 - (1 - i)/2 + 3 * eps, 
		          y = yoff + (1 - i) * height + 5 * eps, 
			  label = i * scale, 
			  rot = -120,
                          gp = gpar(col = labels_color),
		       default.units="native")
                grid.text(x = xoff + i + eps, 
		          y = yoff + -0.05, 
			  label = (1 - i) * scale, 
			  vjust = 1, 
			  rot = 120, 
			  gp = gpar(col = labels_color),
		       default.units="native")
            }
        }


   
 
    xp <- rep(NA,nrow(x))
    yp <- rep(NA,nrow(x))
    for (i in 1:nrow(x))
    {
       xp[i] = as.double(M32[1,] %*% x[i,])
       yp[i] = as.double(M32[2,] %*% x[i,])
    }


    size = unit(if (prop_size)
                    prop_size * (s/max(s))
                else cex, "lines")
		
    grid.points(x = xoff + xp, 
                y = yoff + yp, 
		pch = pch, 
		gp = gpar(col = col), 
		default.units = "native",
                size = size, ...)
	
    if (!is.null(id))
        grid.text(x = xoff + xp, 
	          y = yoff + yp - 0.015, 
		  default.units="native",
		  label = as.character(id), 
		  gp = gpar(col = id_color))
    if (pop)
        popViewport(2)
    else upViewport(2)    
}
