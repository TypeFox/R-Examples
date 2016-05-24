library('likeLTD')
library('ggplot2')
library('stats')

# Case we are going to be looking at.
caseName = 'hammer'
datapath = file.path(system.file("extdata", package="likeLTD"), caseName)

args = list(
  databaseFile = NULL,
  cspFile    = file.path(datapath, 'hammer-CSP.csv'),
  refFile      = file.path(datapath, 'hammer-reference.csv'),
  nUnknowns    = 1,
  doDropin     = TRUE,
  ethnic       = "NDU1",
  adj          = 1.0,
  fst          = 0.02,
  relatedness  = c(0, 0)/4
)

# Create hypothesis for defence
defenceHyp = do.call(defence.hypothesis, args)

graph <- plotLikelihood.2d(defenceHyp, which=c(1,3), x=(1:20/10.0),
                           y=(1:20/10.0))                            + 
                  geom_tile(aes(fill=z))                             + 
                  stat_contour()                                     +
                  labs(title="Log-likelihood map")                   +
                  ylab("Relative contribution of \"Victim 2\"")      +
                  xlab("Relative contribution of unprofiled contributor") 
print(graph)
