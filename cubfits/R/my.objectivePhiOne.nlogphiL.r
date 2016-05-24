### - nlogphiL method will return the negative log (phi * L) value for given
###   phi. These are only called by the optim function to find the MLE.
###
### - These functions are based on per gene list in the following order
###     genes -> amino acids -> synonymous codons -> counts
###   for reu13.list.g, y.g, and n.g.


### Return log likelihood of multinomial distribution + log(phi) value.
### For ROC + NSEf model.
my.objectivePhiOne.nlogphiL.rocnsef <- function(phi, fitlist, reu13.list.g,
    y.g, n.g){
  ret <- my.objectivePhiOne.nlogL.rocnsef(phi, fitlist, reu13.list.g, y.g, n.g)
  ret - log(phi)
} # End of my.objectivePhiOne.nlogphiL.rocnsef().

### For ROC model.
my.objectivePhiOne.nlogphiL.roc <- function(phi, fitlist, reu13.list.g, y.g,
    n.g){
  ret <- my.objectivePhiOne.nlogL.roc(phi, fitlist, reu13.list.g, y.g, n.g)
  ret - log(phi)
} # End of my.objectivePhiOne.nlogphiL.roc().

### For NSEf model.
my.objectivePhiOne.nlogphiL.nsef <- function(phi, fitlist, reu13.list.g, y.g,
    n.g){
  ret <- my.objectivePhiOne.nlogL.nsef(phi, fitlist, reu13.list.g, y.g, n.g)
  ret - log(phi)
} # End of my.objectivePhiOne.nlogphiL.nsef().

