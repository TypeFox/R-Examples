####################################################################################################################
## Author: Ole Christian Lingjaerde
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the clusterGenomics package
## Reference: "Identifying clusters in genomics data by recursive partitioning", Nilsen et al. (2013, preprint)
####################################################################################################################


# Function for plotting heatmap
plotHeatmap = function(X, margin, fast=FALSE, colorscale="green-black-red", pch=".", cex=10) {
   if (missing(margin)) {
      delta = 0.05
      delta2 = 0.4
    } else {
      delta = margin[1]
      delta2 = margin[2]
    }
	plot(0,0,xlim=c(0-delta-delta2,1+delta),ylim=c(0-delta,1+delta+delta2),
		type="n",xaxt="n",yaxt="n",xlab="",ylab="")
	# Basic x-y grid
	Nrow = nrow(X)
	Ncol = ncol(X)
	x0 = seq(0,1,length=Ncol)
	y0 = seq(0,1,length=Nrow)
	# Expand to fill Nrow x Ncol grid
	x1 = rep(x0,Nrow)
	y1 = rep(y0,rep(Ncol,Nrow))
	z1 = as.numeric(t(X))
	# Transform values to colors 
	z1.neg = ifelse(z1<0,-z1,0)
	z1.pos = ifelse(z1>0, z1,0)
	z1.neg.norm = z1.neg / max(abs(z1))
	z1.pos.norm = z1.pos / max(abs(z1))
	if (colorscale == "blue-white-red") {
	  c1 = rgb(1-z1.neg.norm, 1-z1.neg.norm-z1.pos.norm, 1-z1.pos.norm)
	} else if (colorscale == "green-white-red") {
	  c1 = rgb(1-z1.neg.norm, 1-z1.pos.norm, 1-z1.neg.norm-z1.pos.norm)
	} else if (colorscale == "red-white-green") {
	  c1 = rgb(1-z1.pos.norm, 1-z1.neg.norm, 1-z1.neg.norm-z1.pos.norm)
	} else if (colorscale == "green-black-red") {
	  c1 = rgb(z1.pos.norm, z1.neg.norm, 0)
	} else if (colorscale == "red-black-green") {
	  c1 = rgb(z1.neg.norm, z1.pos.norm, 0)
	}
	# Plot heatmap
	if (fast) {
	  plot.xy(xy.coords(x1,y1),type="p",pch=pch,cex=cex, col=c1)
	} else {
	  for (i in 1:Nrow) {
	    for (j in 1:Ncol) {
	      polygon(c(j-1,j-1,j,j)/Ncol, c(i-1,i,i,i-1)/Nrow, col=c1[(i-1)*Ncol+j], border=FALSE)
	    }
	  }
	}
	# Return colorscale
	invisible(list(colorscale=colorscale, value.range=range(z1), 
	color.range=c(max(z1.neg.norm),max(z1.pos.norm))))
}