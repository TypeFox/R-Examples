####################################################################################################################
## Author: Ole Christian Lingjaerde
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the clusterGenomics package
## Reference: "Identifying clusters in genomics data by recursive partitioning", Nilsen et al. (2013, preprint)
####################################################################################################################

# Function for plotting column dendrogram on top
plotTreeCol = function(clust, groups, margin) {
  if (missing(groups)) {
    groups = rep(0,length(clust$order))
  }
  if (missing(margin)) {
    delta = 0.05
    delta2 = 0.4
  } else {
    delta = margin[1]
    delta2 = margin[2]
  }
  invisible(plotForkCol(length(clust$height), clust, as.numeric(groups), delta, delta2))
}



plotForkCol = function(i, clust, groups, delta, delta2) {
  N = length(clust$height)
  mrg = clust$merge
  hgt = clust$height[i]/max(clust$height)
  a = rep(0,2)
  group = rep(0,2)
  b1 = 1 + delta + delta2*hgt    
  for (j in 1:2) {
    if (mrg[i,j]<0) {
      k = which(clust$order == -mrg[i,j])
      a[j] = (k-1)/N
      b0 = 1+delta
      group[j] = groups[k]
    } else {
      tmp = plotForkCol(mrg[i,j],clust,groups,delta,delta2)
      a[j] = tmp[1]
      b0 = 1 + delta + delta2*tmp[2]
      group[j] = tmp[3]
    }
    if (group[j] != -1) {
      lines(c(a[j],a[j]),c(b0,b1),col=group[j]+1)
    } else {
      lines(c(a[j],a[j]),c(b0,b1))
    }
  }
  if (diff(group)==0 & group[1] != -1) {
    lines(c(a[1],a[2]),c(b1,b1),col=group[1]+1)
    invisible(return(c(mean(a),hgt,group[1])))
  } else {
    lines(c(a[1],a[2]),c(b1,b1))
    invisible(return(c(mean(a),hgt,-1)))  
  }
}
