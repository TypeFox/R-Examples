"plotmvoutlier" <-
function (coord, data, quan = 1/2, alpha = 0.025, symb = FALSE, 
    bw=FALSE, plotmap = TRUE, map = "kola.background", which.map = c(1, 
        2, 3, 4), map.col = c(5, 1, 3, 4), map.lwd = c(2, 1, 
        2, 1), pch2=c(3,21), cex2=c(0.7,0.2),col2=c(1,1), 
        lcex.fac=1, ...) 
{
# Mutlivariate outlier plot:
#
# pch2, cex2, col2 ... definitions of symbold for symb=FALSE
# bw ... if TRUE, symbols are in gray-scale (only if symb=TRUE)
# lcex.fac ... factor for multiplication of symbol size (only if symb=TRUE)

    if (ncol(coord) != 2) 
        stop("argument coord has to be two-dimensional")
    rob <- covMcd(data, alpha = quan)
    dist <- mahalanobis(data, center = rob$center, cov = rob$cov)
    xarw <- arw(data, rob$center, rob$cov, alpha = alpha)
    if (xarw$cn != Inf) {
        alpha <- sqrt(c(xarw$cn, qchisq(c(0.75, 0.5, 0.25), ncol(data))))
    }
    else {
        alpha <- sqrt(qchisq(c(0.975, 0.75, 0.5, 0.25), ncol(data)))
    }

    # identify outliers
    rd <- sqrt(dist)
    o <- (rd > sqrt(xarw$cn))

    if (symb == FALSE) {
        if (plotmap == TRUE){ 
            plot(coord, type="n", ...)
            plotbg(map = map, which.map = which.map, map.col = map.col, 
                map.lwd = map.lwd, add.plot = TRUE)
            points(coord[o,], pch=pch2[1],cex=cex2[1],col=col2[1])
            points(coord[!o,], pch=pch2[2],cex=cex2[2],col=col2[2])
	}
        l <- list(outliers = o, md = rd)
    }
    if (symb == TRUE) {
        lpch <- c(3, 3, 16, 1, 1)
        lcex <- c(1.5, 1, 0.5, 1, 1.5)*lcex.fac
        lalpha <- length(alpha)
        xs <- scale(data) - min(scale(data))
        eucl <- sqrt(apply(xs^2, 1, sum))
        if (bw==TRUE){
          rbcol <- rev(gray(seq(from=0,to=0.8,length=nrow(data))))[as.integer(cut(eucl,
            nrow(data), labels = 1:nrow(data)))]
        }
        else{
        rbcol <- rev(rainbow(nrow(data), start = 0, end = 0.7))[as.integer(cut(eucl, 
            nrow(data), labels = 1:nrow(data)))]
        }
	if (plotmap == TRUE){ 
        for (j in 1:lalpha) {
            if (j == 1) {
                plot(coord, type = "n", ...)
        	    plotbg(map = map, which.map = which.map, map.col = map.col, 
             		   map.lwd = map.lwd, add.plot = TRUE)
              	    points(coord[rd >= alpha[j], ], pch = lpch[j], 
                 	   cex = lcex[j], col = rbcol[rd >= alpha[j]])
            }
            if (j > 1 & j < lalpha) 
                points(coord[rd < alpha[j - 1] & rd >= alpha[j], 
                  ], cex = lcex[j], pch = lpch[j], col = rbcol[rd < 
                  alpha[j - 1] & rd >= alpha[j]])
            if (j == lalpha) {
                points(coord[rd < alpha[j - 1] & rd >= alpha[j], 
                  ], cex = lcex[j], pch = lpch[j], col = rbcol[rd < 
                  alpha[j - 1] & rd >= alpha[j]])
                points(coord[rd < alpha[j], ], pch = lpch[j + 
                  1], cex = lcex[j + 1], col = rbcol[rd < alpha[j]])
            }
	}
        }
        l <- list(outliers = o, md = rd, euclidean = eucl)
    }
    l
}
