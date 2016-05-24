####################################################################################################################
## Author: Ole Christian Lingjaerde
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the clusterGenomics package
## Reference: "Identifying clusters in genomics data by recursive partitioning", Nilsen et al. (2013, preprint)
####################################################################################################################

# Function for plotting color bar representing row groups
plotColorbarRow = function(groups, margin) {
  if (missing(margin)) {
    delta = 0.05
    delta2 = 0.4
  } else {
    delta = margin[1]
    delta2 = margin[2]
  }
  groups = as.numeric(groups)
  ugroups = unique(groups)
  for (i in 1:length(groups)) {
    coli = which(groups[i]==ugroups)
    a0 = c(-0.1*delta,-0.9*delta,-0.9*delta,-0.1*delta)
    b0 = c(i-1,i-1,i,i)/length(groups)
    polygon(a0,b0,col=coli,border=FALSE)
  }
}