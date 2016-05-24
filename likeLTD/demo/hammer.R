library('stats')
library('likeLTD')


# Case we are going to be looking at.
caseName = 'hammer'
datapath = file.path(system.file("extdata", package="likeLTD"), caseName)

args = list(
  databaseFile = NULL,
  cspFile    = file.path(datapath, 'hammer-CSP.csv'),
  refFile      = file.path(datapath, 'hammer-reference.csv'),
  nUnknowns    = 0,
  doDropin     = TRUE,
  ethnic       = "NDU1",
  adj          = 1.0,
  fst          = 0.02,
  relatedness  = c(0, 0)/4
)

# Create hypothesis for defence and prosecution.
prosecutionHyp = do.call(prosecution.hypothesis, args)
defenceHyp     = do.call(defence.hypothesis, args)

# Create optimisation parameters, with some modification to hypothesis
# arguments. nUnknowns is modified here and now. It does not need to be, but it
# can be. The somethingParams values contain the objective functions and
# optimisation parameters.
prosecutionParams <- optimisation.params(prosecutionHyp, verbose=TRUE,
                                         nUnknowns=1) 
defenceParams <- optimisation.params(defenceHyp, verbose=FALSE, nUnknowns=2)


# Now perform actual optimization.
prosecutionResult <- do.call(optim, prosecutionParams)
print(prosecutionResult)
# defenceResult     <- do.call(optim, defenceParams)
