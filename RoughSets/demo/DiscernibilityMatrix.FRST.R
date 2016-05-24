########################################################
## Example 1: decision table containing nominal values
########################################################
 library(RoughSets)
 
 data(RoughSetData)
 decision.table <- RoughSetData$hiring.dt 
 
 control.1 <- list(type.relation = c("tolerance", "eq.1"), 
                 type.aggregation = c("t.tnorm", "min"), 
                 t.implicator = "lukasiewicz", type.LU = "implicator.tnorm")
 res.1 <- BC.discernibility.mat.FRST(decision.table, type.discernibility = "standard.red", 
                                     control = control.1)
 
 control.2 <- list(epsilon = 0, delta = 2)
 res.2 <- BC.discernibility.mat.FRST(decision.table, type.discernibility = "gaussian.red",
                                     control = control.2)

 control.3 <- list(type.relation = c("tolerance", "eq.1"), 
                 type.aggregation = c("t.tnorm", "min"),
                 t.implicator = "lukasiewicz", alpha.precision = 0.05)
 res.3 <- BC.discernibility.mat.FRST(decision.table, type.discernibility = "alpha.red", 
                                     control = control.3)

 control.4 <- list(type.relation = c("tolerance", "eq.1"), 
                 type.aggregation = c("t.tnorm", "lukasiewicz"),
                 t.implicator = "lukasiewicz", type.LU = "implicator.tnorm")
 res.4 <- BC.discernibility.mat.FRST(decision.table, type.discernibility = "min.element", 
                                     control = control.4)