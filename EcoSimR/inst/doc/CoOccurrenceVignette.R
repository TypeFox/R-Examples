## ----, warning=FALSE, message=FALSE--------------------------------------
library(EcoSimR)    # load EcoSimR library
set.seed(56)        # for reproducible results

## ----, echo=FALSE, results='asis'----------------------------------------
knitr::kable(dataWiFinches[,1:6], caption='Occurrence matrix of West Indies finches. Data from Gotelli & Abele (1982)')

## ----, echo=FALSE, results='asis'----------------------------------------
Algorithm <- (paste("sim",1:10,sep=""))
E <- "Equiprobable"
P <- "Proportional"
f <- "Fixed"

RowSums <- c(E,f,E,f,P,E,P,P,f,"External Weights")
ColSums <- c(E,E,f,P,f,P,E,P,f, "External Weights")
Notes <- rep("not recommended",9)
Notes[2] <- "recommended"
Notes[9] <- "recommended (DEFAULT)"
Notes[10] <- "recommended"
AlgoTable <- as.data.frame(cbind(Algorithm,
             RowSums, ColSums, Notes))
knitr::kable(AlgoTable, caption='Randomization algorithms for co-occurrence (Gotelli 2000).')


## ----, echo=FALSE,fig.height=4,fig.width=4,fig.align='center'------------
myModel <- cooc_null_model(dataWiFinches,suppressProg=TRUE)
plot(myModel,type="hist")


## ----, echo=FALSE,message=FALSE,warning=FALSE,fig.height=4,fig.width=6,fig.align='center'----
myModel <- cooc_null_model(dataWiFinches,suppressProg=TRUE)
plot(myModel,type="cooc")


## ----, echo=FALSE,message=FALSE,warning=FALSE,fig.height=4,fig.width=4,fig.align='center'----
plot(myModel,type='burn_in')

## ----, eval=FALSE--------------------------------------------------------
#  
#  speciesData          # user must supply a data frame; speciesData=dataWiFinches for default run
#  algo = "sim9"        # randomize occurrences but preserve observed row and columns sums
#  metric = "c_score"   # Stone and Roberts (1990) C-score
#  nReps = 1000         # number of null assemblage created
#  rowNames=TRUE        # reads speciesData as a data frame wtih row labels in the first column
#  saveSeed=FALSE       # if TRUE, saves random number seed
#  burn_in=500          # number of burn-in iterations for sim9
#  algoOpts=list()      # list of other specific options for the algorithm (used for sim10)
#  metricOpts=list()    # list of other specific options for the metric
#  suppressProg= FALSE  # suppress printing of progress bar (for creating markdown files)

## ----,fig.height=4,fig.width=4,fig.align='center'------------------------
# run default settings and show all output
myModel <- cooc_null_model(speciesData=dataWiFinches,suppressProg=TRUE)
summary(myModel)
plot(myModel,type = "hist")

## ----,fig.height=4,fig.width=6,fig.align='center'------------------------
plot(myModel,type = "cooc") 

## ----,fig.height=4,fig.width=4,fig.align='center'------------------------
plot(myModel,type = "burn_in")

# create a model with sim10 and user-supplied species and site weights
myModel <- cooc_null_model(speciesData=dataWiFinches,algo="sim10",
                           suppressProg=TRUE,algoOpts=list(rowWeights
                           =(1:17),colWeights=(1:19)))
summary(myModel)
plot(myModel,type="hist")


## ----,fig.height=4,fig.width=6,fig.align='center'------------------------
plot(myModel,type="cooc")

