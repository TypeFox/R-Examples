 ##############################################
 ## Example: Regression Task
 ##############################################
 library(RoughSets)
 
 data(RoughSetData)
 decision.table <- RoughSetData$housing7.dt 							 
							 
 ## using RI.hybrid.FRST for generating rules
 control <- list(type.aggregation = c("t.tnorm", "lukasiewicz"), t.similarity ="eq.1", 
                  t.implicator = "lukasiewicz")
 rules <- RI.hybridFS.FRST(decision.table, control)
 summary(rules)
 
 ## in this case, we are using the same dataset as training data
 res.1 <- predict(rules, decision.table[, -ncol(decision.table)])
 print("Using RI.hybrid.FRST:")
 print(res.1)
 
 summary(rules)