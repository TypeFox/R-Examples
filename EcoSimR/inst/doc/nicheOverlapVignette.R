## ----,echo=FALSE,warning=FALSE,message=FALSE-----------------------------
library(EcoSimR)

## ----, echo=FALSE, results='asis'----------------------------------------
knitr::kable(dataMacWarb,caption="MacArthur's (1958) warbler data.")

## ----, fig.show='hold', fig.align='center',fig.height=4,fig.width=4,echo=FALSE----
set.seed(56)                              # for repeatable results
myModel <- niche_null_model(dataMacWarb,suppressProg=TRUE)  # default model settings
plot(myModel,type="hist")

## ----,fig.show='hold',fig.height=6,fig.width = 4,fig.align='center'------
plot(myModel,type="niche")

## ----, fig.align='center',echo=FALSE, eval=FALSE-------------------------
#  set.seed(56)                              # for repeatable results
#  myModel <- niche_null_model(dataMacWarb,suppressProg=TRUE)  # default model settings
#  #plot(myModel,type="hist")  #<- throws error, figure margins too large
#  
#  
#  

## ----, eval=FALSE--------------------------------------------------------
#  
#  speciesData          # user must supply a data frame; speciesData=dataMacWarb for default run
#  algo = "ra3"           # reshuffle elements within each row of the matrix
#  metric = "pianka"      # pianka niche overlap index
#  nReps = 1000         # number of null assemblage created
#  rowNames=TRUE        # reads speciesData as a data frame wtih row labels in the first column
#  saveSeed=FALSE       # if TRUE, saves random number seed
#  algoOpts=list()      # list of other specific options for the algorithm
#  metricOpts=list()    # list of other specific options for the metric
#  suppressProg= FALSE  # suppress printing of progress bar (for creating markdown files)

## ----, eval=FALSE--------------------------------------------------------
#  str(dataMacWarb)  # structure of MacArthur's warbler data set
#  summary(myModel)  # output summary of null model analysis
#  
#  #create a random data set with uniform (0,1) values
#  myRandomData <- matrix(runif(300), nrow=30)
#  
#  # run null model with czekanowski index and ra1, 5000 replications
#  myRandomModel <- niche_null_model(speciesData=myRandomData, rowNames=FALSE,
#                              algo="ra1", metric="czekanowski",
#                              suppressProg=TRUE,nReps=5000)
#  
#  # print summary of model and plot histogram
#  summary(myRandomModel)
#  plot(myRandomModel,type="hist")
#  

