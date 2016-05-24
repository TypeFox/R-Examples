 #################################################################
 ## Example: Determine cut values and generate new decision table
 #################################################################
 library(RoughSets)
 
 dt.ex1 <- data.frame(c(1, 1.2, 1.3, 1.4, 1.4, 1.6, 1.3), c(2, 0.5, 3, 1, 2, 3, 1),
                              c(1, 0, 0, 1, 0, 1, 1))
 colnames(dt.ex1) <- c("a", "b", "d")
 decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 3, 
                                      indx.nominal = c(3)) 

 cut.values <- D.discretize.equal.intervals.RST(decision.table, nOfIntervals = 4)

 ## generate new decision table
 new.decTable <- SF.applyDecTable(decision.table, cut.values)