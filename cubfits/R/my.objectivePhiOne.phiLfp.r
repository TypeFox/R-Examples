### - PM method will return the phi * L * 1 for given phi
###   These are only called by the integrate function to find the posterior
###   mean.
### - These functions are based on per gene list in the following order
###     genes -> amino acids -> synonymous codons -> counts
###   for reu13.list.g, y.g, and n.g.
### - add.llv.g is to scale log likelihood value for a gene.


### Return x * Likelihood(phi) * flat prior
### For ROC + NSEf model.
my.objectivePhiOne.xLfp.rocnsef <- function(phi, fitlist, reu13.list.g, y.g,
    n.g, add.llv.g = 0.0){
  ret <- my.objectivePhiOne.nlogL.rocnsef(phi, fitlist, reu13.list.g, y.g, n.g)
  ret <- exp(-ret + log(phi) + add.llv.g)
  ret
} # End of my.objectivePhiOne.xLfp.rocnsef().

### For ROC model.
my.objectivePhiOne.xLfp.roc <- function(phi, fitlist, reu13.list.g, y.g, n.g,
    add.llv.g = 0.0){
  ret <- my.objectivePhiOne.nlogL.roc(phi, fitlist, reu13.list.g, y.g, n.g)
  ret <- exp(-ret + log(phi) + add.llv.g)
  ret
} # End of my.objectivePhiOne.xLfp.roc().

### For NSEf model.
my.objectivePhiOne.xLfp.nsef <- function(phi, fitlist, reu13.list.g, y.g, n.g,
    add.llv.g = 0.0){
  ret <- my.objectivePhiOne.nlogL.nsef(phi, fitlist, reu13.list.g, y.g, n.g)
  ret <- exp(-ret + log(phi) + add.llv.g)
  ret
} # End of my.objectivePhiOne.xLfp.nsef().

