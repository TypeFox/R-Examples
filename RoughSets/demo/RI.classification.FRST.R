 ##############################################
 ## Example: Classification Task
 ##############################################
 library(RoughSets)
 
 data(RoughSetData)
 decision.table <- RoughSetData$pima7.dt 							 
							 
 ## using RI.hybrid.FRST for generating rules
 control <- list(type.aggregation = c("t.tnorm", "lukasiewicz"), t.similarity ="eq.3", 
                 t.implicator = "lukasiewicz")
 rules.1 <- RI.hybridFS.FRST(decision.table, control)

 ## in this case, we are using the same dataset as training data
 res.1 <- predict(rules.1, decision.table[, -ncol(decision.table), drop = FALSE])

 ## using RI.GFRS.FRST for generating rules
 control <- list(alpha.precision = 0.05, type.aggregation = c("t.tnorm", "lukasiewicz"), 
                 t.similarity ="eq.3", t.implicator = "lukasiewicz")						 
 rules.2 <- RI.GFRS.FRST(decision.table, control)

 ## in this case, we are using the same dataset as training data
 res.2 <- predict(rules.2, decision.table[, -ncol(decision.table)])

 print("Results:")
 print("Using RI.hybrid.FRST:")
 print(res.1)
 print("Using RI.GFRS.FRST:")
 print(res.2)
 
 summary(rules.1)
 summary(rules.2)