plot.fds = function (x, plot.type = c("functions", "time", "depth", "density"), 
    col = NULL, type = "l", lty = 1, xlab = x$xname, ylab = x$yname, 
    pch = c(1:9, 0, letters, LETTERS), add = FALSE, index, 
    colorchoice = c("rainbow", "heat.colors","terrain.colors","topo.colors","cm.colors",
					"rainbow_hcl", "gray", "sequential_hcl", "heat_hcl",
					"terrain_hcl", "diverge_hcl"), plotlegend = FALSE, legendpos = "topright",
					ncol = 1, ...) 
{
    colorchoice = match.arg(colorchoice)
    if (class(x)[1] == "fts" | class(x)[1] == "fds" | class(x)[1] == "sfts") 
    {
        plot.type <- match.arg(plot.type)
        if (plot.type == "time") {
            if (class(x)[1] == "fts" | class(x)[1] == "sfts") {
                if (is.null(col)) {
                  nx <- length(x$x)
                  col = rainbow(min(1024, 1.25 * nx))
                }
                if (xlab == x$xname) {
                  xlab <- "Time"
                }
                if(is.null(colnames(x$y)))
                {
                    year = 1:ncol(x$y)
                }
                else
                {
                    year = as.numeric(colnames(x$y))
                }
                if (add == FALSE) {
                  matplot(year, t(x$y), type = type, ylab = ylab, 
                    xlab = xlab, col = col, lty = lty, pch = pch, 
                    ...)
                }
                else {
                  matlines(year, t(x$y)[, index], type = type, 
                    ylab = ylab, xlab = xlab, col = col, lty = lty, 
                    pch = pch, ...)
                }
            }
            else {
                stop("object is not a functional time series.")
            }
        }
        else {  
            if (is.null(col)) {
                ny <- ncol(as.matrix(x$y))
                if (ny > 1) 
				{
					if(colorchoice == "rainbow")
					{
                    	col <- rainbow(min(1024, 1.25 * ny))
					}				
					if(colorchoice == "heat.colors")
					{
						col <- heat.colors(min(1024, 1.25 * ny))
					}
					if(colorchoice == "terrain.colors")
					{
						col <- terrain.colors(min(1024, 1.25 * ny))
					}
					if(colorchoice == "topo.colors")
					{
						col <- topo.colors(min(1024, 1.25 * ny))
					}
					if(colorchoice == "cm.colors")
					{
						col <- cm.colors(min(1024, 1.25 * ny))
					}
					if(colorchoice == "rainbow_hcl")
					{
                		col <- rainbow_hcl(min(1024, 1.25 * ny))
					}
					if(colorchoice == "gray")
					{
                        col <- gray(ny:1/ny)
					}
					if(colorchoice == "sequential_hcl")
					{
						col <- sequential_hcl(min(1024, 1.25 * ny))
					}
					if(colorchoice == "heat_hcl")
					{
						col <- heat_hcl(min(1024, 1.25 * ny))
					}
					if(colorchoice == "terrain_hcl")
					{
						col <- terrain_hcl(min(1024, 1.25 * ny))
					}
					if(colorchoice == "diverge_hcl")
					{
						col <- diverge_hcl(min(1024, 1.25 * ny))
					}	
				}
				else {
                  col <- 1
                }
            }
            yy <- as.matrix(x$y)
            if (plot.type == "depth") {
                sco <- PCAproj(t(yy), k = 2, center = median)$score
                center <- compute.bagplot(sco)$center
                lineindex <- order(mahalanobis(sco, center, cov(sco)))
                yy <- yy[, lineindex]
                yymax <- yy[, 1]
            }
            else if (plot.type == "density") {
                sco <- PCAproj(t(yy), k = 2, center = median)$score
                X <- cbind(sco[, 1], sco[, 2])
                h = Hscv.diag(X, binned = TRUE)
                den = kde(x = X, H = h)
                den = list(x = den$eval.points[[1]], y = den$eval.points[[2]], 
                  z = den$estimate)
                den2 <- hdrcde::hdr.2d(sco[, 1], sco[, 2], prob = c(0.01, 0.5), den)
                lineindex <- order(den2$fxy, decreasing = TRUE)
                yy <- yy[, lineindex]
                yymax <- yy[, 1]
            }
            if (nrow(yy) == 1) 
                plot(ts(c(yy), start = start(x$time), frequency = frequency(x$time)), 
                  ylab = ylab, ...)
            if (plot.type == "functions") {
                if (add == FALSE) {
                  matplot(x$x, yy, col = col, xlab = xlab, ylab = ylab, 
                    type = type, lty = lty, pch = pch, ...)
					if(plotlegend == TRUE)
					{
						legend(legendpos,c(colnames(x$y)[1], 
								colnames(x$y)[ceiling(ny/2)],
								colnames(x$y)[length(colnames(x$y))]),
								col=c(col[1],col[ceiling(ny/2)],col[ny]), ncol=ncol, lty=1)			
					}
                }
                else {
                  matlines(x$x, yy[, index], col = col, xlab = xlab, 
                    ylab = ylab, type = type, lty = lty, pch = pch, 
                    ...)
                }
            }
            else {
					if (add == FALSE) 
					{
						matplot(x$x, yy, col = col, xlab = xlab, ylab = ylab, 
								type = type, lty = lty, pch = pch, ...)
						lines(x$x, yymax, col = "black", ...)
						if(plotlegend == TRUE)
						{
							legend(legendpos,c(colnames(yy)[1], 
									colnames(yy)[ceiling(ny/2)],
									colnames(yy)[length(colnames(yy))]),
									col=c(col[1],col[ceiling(ny/2)],col[ny]), ncol=ncol,lty=1)		
						}
					}
                else {
                  matlines(x$x, yy[, index], col = col, xlab = xlab, 
                    ylab = ylab, type = type, lty = lty, pch = pch, 
                    ...)
                }
            }
        }
    }
    else {
        stop("object is not a functional model.")
    }
}
