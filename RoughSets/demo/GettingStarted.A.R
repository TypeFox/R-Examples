##############################################################
## A.1 Example: Basic concepts of rough set theory
##############################################################
## Using hiring data set, see RoughSetData
data(RoughSetData)
decision.table <- RoughSetData$hiring.dt

## define considered attributes which are first, second, and
## third attributes
attr.P <- c(1,2,3)

## compute indiscernibility relation
IND <- BC.IND.relation.RST(decision.table, feature.set = attr.P)

## compute lower and upper approximations
roughset <- BC.LU.approximation.RST(decision.table, IND)

## Determine regions
region.RST <- BC.positive.reg.RST(decision.table, roughset)

## The decision-relative discernibility matrix and reduct
disc.mat <- BC.discernibility.mat.RST(decision.table, range.object = NULL)

##############################################################
## A.2 Example: Basic concepts of fuzzy rough set theory
##############################################################
## Using pima7 data set, see RoughSetData
data(RoughSetData)
decision.table <- RoughSetData$pima7.dt

## In this case, let us consider the first and second attributes
conditional.attr <- c(1, 2)

## We are using the "lukasiewicz" t-norm and the "tolerance" relation
## with "eq.1" as fuzzy similarity equation
control.ind <- list(type.aggregation = c("t.tnorm", "lukasiewicz"),
type.relation = c("tolerance", "eq.1"))

## Compute fuzzy indiscernibility relation
IND.condAttr <- BC.IND.relation.FRST(decision.table, attributes = conditional.attr,
control = control.ind)

## Compute fuzzy lower and upper approximation using type.LU : "implicator.tnorm"
## Define index of decision attribute
decision.attr = c(9)

## Compute fuzzy indiscernibility relation of decision attribute
## We are using "crisp" for type of aggregation and type of relation
control.dec <- list(type.aggregation = c("crisp"), type.relation = "crisp")
IND.decAttr <- BC.IND.relation.FRST(decision.table, attributes = decision.attr,
control = control.dec)

## Define control parameter containing type of implicator and t-norm
control <- list(t.implicator = "lukasiewicz", t.tnorm = "lukasiewicz")

## Compute fuzzy lower and upper approximation
FRST.LU <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr,
type.LU = "implicator.tnorm", control = control)

## Determine fuzzy positive region and its degree of dependency
fuzzy.region <- BC.positive.reg.FRST(decision.table, FRST.LU)

###############################################################
## B Example : Data analysis based on RST and FRST
## In this example, we are using wine dataset for both RST and FRST
###############################################################
## Using wine data set, see RoughSetData
## In this case, we use data at indices 6 - 178 as the training data
## Then, we define the column 14 as the decision and nominal attribute
## Not run: data(RoughSetData)
wine.decTable <- SF.asDecisionTable(dataset = RoughSetData$wine.dt[6 : 178, ],
decision.attr = 14, indx.nominal = 14)

## define first five instances as newdata or testing data
tst.wine <- SF.asDecisionTable(dataset = wine.decTable[1:5, -ncol(wine.decTable)])

###############################################################
## B.1 Example : Rough Set Theory
###############################################################
## DISCRETIZATION STEP
## In this example, we are using local strategy algorithm
cut.values.tra <- D.discretization.RST(wine.decTable,
type.method = "global.discernibility")

## generate new decision table
d.new.tra.rst <- SF.applyDecTable(wine.decTable, cut.values.tra)
d.new.tst.rst <- SF.applyDecTable(tst.wine, cut.values.tra)

## FEATURE SELECTION STEP
## In this example, we are using permutation algorithm
## which generates a single superreduct
red.rst <- FS.feature.subset.computation(d.new.tra.rst,
method = "quickreduct.rst")

## generate new decision table according to the reduct (optional)
fs.new.decTable.rst <- SF.applyDecTable(d.new.tra.rst, red.rst)

## RULE INDUCTION
## In this case, we are using the original decision table,
## we also can use the decision table resulting from feature selection
rules.rst <- RI.indiscernibilityBasedRules.RST(d.new.tra.rst, red.rst)

## predicting newdata
res.1 <- predict(rules.rst, d.new.tst.rst)