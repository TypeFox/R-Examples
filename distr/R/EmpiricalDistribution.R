###############################################################################
## Generating function to generate empirical distribution given some data
###############################################################################

## simple wrapper to DiscreteDistribution
EmpiricalDistribution <- function(data, .withArith = FALSE,
                                  .withSim = FALSE, .lowerExact = TRUE, .logExact = FALSE,
                                  .DistrCollapse = 
                                    getdistrOption("DistrCollapse"),
                                  .DistrCollapse.Unique.Warn = 
                                    getdistrOption("DistrCollapse.Unique.Warn"),
                                  .DistrResolution = getdistrOption("DistrResolution"),
                                  Symmetry = NoSymmetry()){
  DiscreteDistribution(supp = data, .withArith = .withArith, .withSim = .withSim,
                       .lowerExact = .lowerExact, .logExact = .logExact, 
                       .DistrCollapse = .DistrCollapse, 
                       .DistrCollapse.Unique.Warn = .DistrCollapse.Unique.Warn,
                       .DistrResolution = .DistrResolution,
                       Symmetry = Symmetry)
}
