## ----eval=FALSE----------------------------------------------------------
#  install.packages("SDR")

## ----eval=FALSE----------------------------------------------------------
#  devtools::install_github('aklxao2/SDR')

## ----eval=FALSE----------------------------------------------------------
#  irisTraining <- read.keel("irisTra.dat")
#  irisTest <- read.keel("irisTst.dat")

## ----eval=FALSE----------------------------------------------------------
#  irisTraining <- read.keel("irisTra.dat", nLabels = 5)
#  irisTest <- read.keel("irisTst.dat", nLabels = 5)

## ----eval=FALSE----------------------------------------------------------
#  irisTraining <- keelFromDataFrame(data = iris, relation = "iris",)

## ----eval=FALSE----------------------------------------------------------
#  irisTraining <- keelFromDataFrame(data = iris, relation = "iris", nLabels = 5, names = c("Sepal.Length, Sepal.Width", "Petal.Length", "Petal.Width"), types = c('r','r','r','r','c'), classNames = c("setosa", "virginica","versicolor"))

## ----eval=FALSE----------------------------------------------------------
#  MESDIF(paramFile = "MESDIFparameters.txt")

## ----highlight=TRUE------------------------------------------------------
library("SDR")
MESDIF( paramFile = NULL,
        training = habermanTra,
        test = habermanTst,
        output = c("optionsFile.txt", "rulesFile.txt", "testQM.txt"),
        seed = 0,
        nLabels = 3,
        nEval = 300,
        popLength = 100,
        eliteLength = 3,
        crossProb = 0.6,
        mutProb = 0.01,
        RulesRep = "can",
        Obj1 = "CSUP",
        Obj2 = "CCNF",
        Obj3 = "null",
        Obj4 = "null",
        targetClass = "positive"
        )

## ----eval=FALSE----------------------------------------------------------
#  SDR_GUI()

