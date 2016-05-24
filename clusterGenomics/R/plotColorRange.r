####################################################################################################################
## Author: Ole Christian Lingjaerde
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the clusterGenomics package
## Reference: "Identifying clusters in genomics data by recursive partitioning", Nilsen et al. (2013, preprint)
####################################################################################################################

# Function for plotting color range
plotColorRange = function(colrange) {
  plot(0,0,type="n",xlim=colrange$value.range,ylim=c(0,1),xlab="",ylab="",xaxt="n",yaxt="n")
  abline(v=0)
  axis(1, colrange$value.range, round(colrange$value.range*100)/100)
  axis(1, 0, 0)
  colvec = rev(seq(0, colrange$color.range[1], length=100))
  posvec = -rev(seq(0, -colrange$value.range[1], length=100))
  if (colrange$colorscale == "blue-white-red") {
	  c1 = rgb(1-colvec, 1-colvec, 1)
  } else if (colrange$colorscale == "green-white-red") {
	  c1 = rgb(1-colvec, 1, 1-colvec)
  } else if (colrange$colorscale == "red-white-green") {
	  c1 = rgb(1, 1-colvec, 1-colvec)
  } else if (colrange$colorscale == "green-black-red") {
	  c1 = rgb(0, colvec, 0)
  } else if (colrange$colorscale == "red-black-green") {
	  c1 = rgb(colvec, 0, 0)
  }  
  for (i in 1:(length(colvec)-1)) {
    polygon(c(posvec[i],posvec[i],posvec[i+1],posvec[i+1]), c(0,1,1,0), 
      col=c1[i], border=FALSE)
  }
  
  colvec = seq(0, colrange$color.range[2], length=100)
  posvec = seq(0, colrange$value.range[2], length=100)
  if (colrange$colorscale == "blue-white-red") {
	  c1 = rgb(1, 1-colvec, 1-colvec)
  } else if (colrange$colorscale == "green-white-red") {
	  c1 = rgb(1, 1-colvec, 1-colvec)
  } else if (colrange$colorscale == "red-white-green") {
	  c1 = rgb(1-colvec, 1, 1-colvec)
  } else if (colrange$colorscale == "green-black-red") {
	  c1 = rgb(colvec, 0, 0)
  } else if (colrange$colorscale == "red-black-green") {
	  c1 = rgb(0, colvec, 0)
  }  
  for (i in 1:(length(colvec)-1)) {
    polygon(c(posvec[i],posvec[i],posvec[i+1],posvec[i+1]), c(0,1,1,0), 
      col=c1[i], border=FALSE)
  }  
  
}
