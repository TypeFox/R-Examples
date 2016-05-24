### - Thses will return the L * 1 for given phi
###   These are only called by the integrate function to find the posterior
###   mean.
### - These functions are based on per gene list in the following order
###     genes -> amino acids -> synonymous codons -> counts
###   for reu13.list.g, y.g, and n.g.
### - add.llv.g is to scale log likelihood value for a gene.


### Return Likelihood(phi) * flat prior
### For ROC + NSEf model.
my.objectivePhiOne.Lfp.rocnsef <- function(phi, fitlist, reu13.list.g, y.g, n.g,
    add.llv.g = 0.0){
  ret <- my.objectivePhiOne.nlogL.rocnsef(phi, fitlist, reu13.list.g, y.g, n.g)
  exp(-ret + add.llv.g)
} # End of my.objectivePhiOne.Lfp.rocnsef().

### For ROC model.
my.objectivePhiOne.Lfp.roc <- function(phi, fitlist, reu13.list.g, y.g, n.g,
    add.llv.g = 0.0){
  ret <- my.objectivePhiOne.nlogL.roc(phi, fitlist, reu13.list.g, y.g, n.g)
  exp(-ret + add.llv.g)
} # End of my.objectivePhiOne.Lfp.roc().

### For NSEf model.
my.objectivePhiOne.Lfp.nsef <- function(phi, fitlist, reu13.list.g, y.g, n.g,
    add.llv.g = 0.0){
  ret <- my.objectivePhiOne.nlogL.nsef(phi, fitlist, reu13.list.g, y.g, n.g)
  exp(-ret + add.llv.g)
} # End of my.objectivePhiOne.Lfp.nsef().

