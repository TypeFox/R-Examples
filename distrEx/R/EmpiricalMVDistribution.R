###############################################################################
## Generating function to generate empirical distribution given some data
###############################################################################

EmpiricalMVDistribution <- function(data, Symmetry = NoSymmetry()){
  DiscreteMVDistribution(supp = data, Symmetry = Symmetry)
}
