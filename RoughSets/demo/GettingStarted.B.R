###############################################################
## B Example : Data analysis based on RST and FRST
## In this example, we are using wine dataset for both RST and FRST
###############################################################
## Using wine data set
data(RoughSetData)
wine.dt <- RoughSetData$wine.dt

## shuffle the data
set.seed(5) 
wine.dt <- wine.dt[sample(nrow(wine.dt)),]

## split the data and then construct the decisionTable format for training
num.percent <- round(0.8*nrow(wine.dt))
wine.tra <- SF.asDecisionTable(dataset = wine.dt[1:num.percent, ],
    decision.attr = 14, indx.nominal = 14)

## check properties of training data
print(attr(wine.tra, "nominal.attrs"))
print(attr(wine.tra, "desc.attrs"))
print(attr(wine.tra, "decision.attr"))

## construct testing data
## we do not define decision.attr
wine.tst <- SF.asDecisionTable(wine.dt[(num.percent+1):
    nrow(wine.dt), -ncol(wine.dt)])

###############################################################
## B.1 Example : Rough Set Theory
###############################################################
## DISCRETIZATION STEP
## In this example, we are using local strategy algorithm
cutValues <- D.discretization.RST(wine.tra,
    type.method = "global.discernibility")
print(cutValues)
	
## generate new decision table
wine.tra.d <- SF.applyDecTable(wine.tra, cutValues)
wine.tst.d <- SF.applyDecTable(wine.tst, cutValues)

## FEATURE SELECTION STEP
## For example, we are using permutation algorithm
## which generates a single superreduct
reduct <- FS.feature.subset.computation(wine.tra.d,
    method = "quickreduct.rst")
print(reduct)

## generate new decision table according to the reduct (optional)
wine.tra.fs <- SF.applyDecTable(wine.tra.d, reduct)

## RULE INDUCTION
## Note: because the method considers reduct as input data, 
## we can supply the original decision table
decRules.rst <- RI.indiscernibilityBasedRules.RST(wine.tra.d, reduct)
summary(decRules.rst)

## prediction
predValues.rst <- predict(decRules.rst, wine.tst.d)

###############################################################
## B.2 Example : Fuzzy Rough Set Theory
###############################################################
## Note: we do not need to conduct discretization

## FEATURE SELECTION
reduct <- FS.feature.subset.computation(wine.tra,
    method = "quickreduct.frst")
print(reduct)

## generate new decision tables
wine.tra.fs <- SF.applyDecTable(wine.tra, reduct)
wine.tst.fs <- SF.applyDecTable(wine.tst, reduct)

## INSTANCE SELECTION
## It should be noted in this case we are using the decision table
## that results from feature selection
indx <- IS.FRIS.FRST(wine.tra.fs, control =
    list(threshold.tau = 0.2, alpha = 1))

## generate a new decision table
wine.tra.is <- SF.applyDecTable(wine.tra.fs, indx)

## CLASSIFIER
## Predict the testing data by using the FRNN.O method
control.frnn.o <- list(m = 2, type.membership = "gradual")

predValues.frnn.o <- C.FRNN.O.FRST(wine.tra.is, newdata = wine.tst.fs,
   control = control.frnn.o)

## FRNN
control.frnn <- list(type.LU = "implicator.tnorm", k = 20, 
                 type.aggregation = c("t.tnorm", "lukasiewicz"), 
                 type.relation = c("tolerance", "eq.1"), t.implicator = "lukasiewicz") 									   
predValues.frnn <- C.FRNN.FRST(wine.tra.is, newdata = wine.tst.fs,
                              control = control.frnn)

## POSNN
control.posnn <- list(type.LU = "implicator.tnorm", k = 20, t.tnorm = "lukasiewicz", 
                 type.relation = c("tolerance", "eq.1"), t.implicator = "lukasiewicz")

predValues.posnn <- C.POSNN.FRST(wine.tra.is, 
                newdata = wine.tst.fs, control = control.posnn)						  
							  
## RULE INDUCTION (Rule-based classifiers)
## In this case, we are using the RI.hybridFS.FRST method
control.ri <- list(type.aggregation = c("t.tnorm", "lukasiewicz"),
type.relation = c("tolerance", "eq.3"), t.implicator = "kleene_dienes")

decRules.hybrid.a <- RI.hybridFS.FRST(wine.tra.is, control.ri)
decRules.hybrid.b <- RI.hybridFS.FRST(wine.tra, control.ri)

## predicting newdata
predValues.hybrid.a <- predict(decRules.hybrid.a, wine.tst.fs)
predValues.hybrid.b <- predict(decRules.hybrid.b, wine.tst)

## check error
realValues <- wine.dt[(num.percent+1):
    nrow(wine.dt), ncol(wine.dt), drop = FALSE]
	
err.rst = 100*sum(realValues != predValues.rst)/nrow(predValues.rst)
err.frnn.o = 100*sum(realValues != predValues.frnn.o)/nrow(predValues.frnn.o)
err.frnn = 100*sum(realValues != predValues.frnn)/nrow(predValues.frnn)
err.posnn = 100*sum(realValues != predValues.posnn)/nrow(predValues.posnn)
err.hybrid.a = 100*sum(realValues != predValues.hybrid.a)/nrow(predValues.hybrid.a)
err.hybrid.b = 100*sum(realValues != predValues.hybrid.b)/nrow(predValues.hybrid.b)

cat("Classification Error of RI based on RST:", err.rst, "\n")
cat("Classification Error of FRNN.0:", err.frnn.o, "\n")
cat("Classification Error of FRNN:", err.frnn)
cat("Classification Error of POSNN:", err.posnn)
cat("Classification Error of RI.hybrid with FS and IS:", err.hybrid.a)
cat("Classification Error of RI.hybrid without FS and IS:", err.hybrid.b)




