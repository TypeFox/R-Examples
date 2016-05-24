## ----, warning=FALSE, message=FALSE--------------------------------------
library(EcoSimR)    # load EcoSimR library
set.seed(56)        # for reproducible results

## ----, echo=FALSE, results='asis'----------------------------------------
knitr::kable(dataRodents, caption='Average body sizes of Sonoran desert rodents. Data from Brown (1975)')

## ----, echo=FALSE, fig.height=4,fig.width=4,fig.align='center'-----------
myModel <- cooc_null_model(dataWiFinches,suppressProg=TRUE)
plot(myModel,type="hist")


## ----, fig.height=8,fig.width=4,fig.align='center'-----------------------
myModel <- size_null_model(dataRodents,suppressProg=TRUE)
plot(myModel,type="size") 

## ----, eval=FALSE--------------------------------------------------------
#  
#  speciesData           # user must supply a data frame; speciesData=dataWiFinches for default run
#  algo = "size_uniform" # randomize interior species with a uniform distribution
#  metric = "var_ratio"  # variance of size ratios of adjacent species
#  nReps = 1000          # number of null assemblage created
#  rowNames=TRUE         # reads speciesData as a data frame wtih row labels in the first column
#  saveSeed=FALSE        # if TRUE, saves random number seed
#  burn_in=500           # number of burn-in iterations for sim9
#  algoOpts=list()       # list of other specific options for the algorithm (used for size_source_pool)
#  metricOpts=list()     # list of other specific options for the metric
#  suppressProg= FALSE   # suppress printing of progress bar (for creating markdown files)

## ----, fig.height=4,fig.width=4,fig.align='center'-----------------------
# run default settings and show all output
myModel <- size_null_model(speciesData=dataRodents,suppressProg=TRUE)
summary(myModel)
plot(myModel,type = "hist")

## ----,fig.height=8,fig.width=4,fig.align='center'------------------------
plot(myModel,type="size") # throws error in vignette: figure margins too large

## ----, fig.height=4,fig.width=4,fig.align='center'-----------------------
# test for minimum size differences with a source pool model

# create a source pool of the rodent body sizes plus 20 other species
mySource <-c(dataRodents$Sonoran,as.double(sample(150,20)))
           
# create an arbitrary set of probabilty weights
myProbs <- runif(26)

# run the model
myModel <- size_null_model(speciesData=dataRodents,suppressProg=TRUE,
           metric="min_diff",algo="size_source_pool",
           algoOpts=list(sourcePool=mySource,speciesProbs=myProbs))

# show the results
summary(myModel)
plot(myModel,type="hist")

## ----, fig.height=8,fig.width=4,fig.align='center'-----------------------
plot(myModel,type="size") 

