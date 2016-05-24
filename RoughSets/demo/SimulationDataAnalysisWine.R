########################################
## Simple Examples
########################################
## The data must be in data.frame
dt.ex1 <- data.frame(
  c(100.2, 102.6, NA, 99.6, 99.8, 96.4, 96.6, NA),
  c(NA, "yes", "no", "yes", NA, "yes", "no", "yes"),
  c("no", "yes", "no", "yes", "yes", "no", "yes", NA),
  c("yes", "yes", "no", "yes", "no", "no", "no", "yes"))

## We should define the names of attributes    
colnames(dt.ex1) <- c("Temp","Headache","Nausea","Flu")

## Construct the DecisionTable format,
## we define that the 2,3,4 attributes are nominal and
## the 4 attribute is the decision attribute
decTable.ex1 <- SF.asDecisionTable(dataset=dt.ex1, 
   decision.att =4, indx.nominal = c(2:4))
print(attr(decTable.ex1, "nominal.attrs"))
print(attr(decTable.ex1, "desc.attrs"))
print(attr(decTable.ex1, "decision.attr"))

## Because we have missing values, the functions of missing value completion should be performed
obj.MV <- MV.missingValueCompletion(decTable.ex1, 
    type.method = "deletionCases")

## Generate a new decision table, which is without missing values     
new.decTable <- SF.applyDecTable(decTable.ex1, obj.MV)
print(new.decTable)
print(attr(new.decTable, "nominal.attrs"))
print(attr(new.decTable, "desc.attrs"))
print(attr(new.decTable, "decision.attr"))

## In order to use functions based on RST we need to perform a discretization of numeric attributes
obj.cutVal <- D.discretization.RST(new.decTable, 
   type.method = "global.discernibility")

## Generate a new decision table using cut values from discretization
new.decTable <- SF.applyDecTable(new.decTable, 
   obj.cutVal)
print(new.decTable) 

#################################################
## Examples: Data analysis using the wine dataset
## 1. Learning and prediction using RST
#################################################
## Load the data
data(RoughSetData)
dataset <- RoughSetData$wine.dt

## Shuffle the data with set.seed
set.seed(5) 
dt.Shuffled <- dataset[sample(nrow(dataset)),]

## Split the data into training and testing
idx <- round(0.8 * nrow(dt.Shuffled))
wine.tra <-SF.asDecisionTable(dt.Shuffled[1:idx,],
 decision.attr = 14, indx.nominal = 14)
wine.tst <- SF.asDecisionTable(dt.Shuffled[
 (idx+1):nrow(dt.Shuffled), -ncol(dt.Shuffled)])
 
## DISCRETIZATION
cut.values <- D.discretization.RST(wine.tra,
 type.method = "global.discernibility")
d.tra <- SF.applyDecTable(wine.tra, cut.values)
d.tst <- SF.applyDecTable(wine.tst, cut.values)

## FEATURE SELECTION
red.rst <- FS.feature.subset.computation(d.tra, 
 method="quickreduct.rst")
fs.tra <- SF.applyDecTable(d.tra, red.rst)
## RULE INDUCTION
rules <- RI.indiscernibilityBasedRules.RST(d.tra, 
 red.rst)

## predicting newdata
pred.vals <- predict(rules, d.tst)

#################################################
## Examples: Data analysis using the wine dataset
## 2. Learning and prediction using FRST
#################################################

## FEATURE SELECTION
reduct <- FS.feature.subset.computation(wine.tra,
 method = "quickreduct.frst")

## generate new decision tables
wine.tra.fs <- SF.applyDecTable(wine.tra, reduct)
wine.tst.fs <- SF.applyDecTable(wine.tst, reduct)

## INSTANCE SELECTION
indx <- IS.FRIS.FRST(wine.tra.fs, 
  control = list(threshold.tau = 0.2, alpha = 1))

## generate a new decision table
wine.tra.is <- SF.applyDecTable(wine.tra.fs, indx)

## RULE INDUCTION (Rule-based classifiers)
control.ri <- list(
 type.aggregation = c("t.tnorm", "lukasiewicz"), 
 type.relation = c("tolerance", "eq.3"), 
 t.implicator = "kleene_dienes")

decRules.hybrid <- RI.hybridFS.FRST(wine.tra.is, 
  control.ri)

## predicting newdata
predValues.hybrid <- predict(decRules.hybrid, 
  wine.tst.fs)
 
#################################################
## Examples: Data analysis using the wine dataset
## 3. Prediction using fuzzy nearest neighbor classifiers
#################################################

## using FRNN.O
control.frnn.o <- list(m = 2, 
  type.membership = "gradual")

predValues.frnn.o <- C.FRNN.O.FRST(wine.tra.is, 
  newdata = wine.tst.fs, control = control.frnn.o)

## Using FRNN
control.frnn <- list(type.LU = "implicator.tnorm",k=20, 
  type.aggregation = c("t.tnorm", "lukasiewicz"), 
  type.relation = c("tolerance", "eq.1"), 
  t.implicator = "lukasiewicz") 
  									   
predValues.frnn <- C.FRNN.FRST(wine.tra.is, 
  newdata = wine.tst.fs, control = control.frnn)

## calculating error
real.val <- dt.Shuffled[(idx+1):nrow(dt.Shuffled),
 ncol(dt.Shuffled), drop = FALSE]
 
err.1 <- 100*sum(pred.vals!=real.val)/nrow(pred.vals)
err.2 <- 100*sum(predValues.hybrid!=real.val)/
  nrow(predValues.hybrid)
err.3 <- 100*sum(predValues.frnn.o!=real.val)/
  nrow(predValues.frnn.o)
err.4 <- 100*sum(predValues.frnn!=real.val)/
  nrow(predValues.frnn)

cat("The error percentage of RI based on RST = ", err.1, "\n")
cat("The error percentage of RI based on FRST = ", err.2, "\n")
cat("The error percentage of FRNN.O = ", err.3, "\n")
cat("The error percentage of FRNN = ", err.4, "\n")
