 ##############################################
 ## Example: Classification Task
 ##############################################
 library(RoughSets)
 
 data(RoughSetData)
 decision.table <- RoughSetData$hiring.dt 							 
		
 #### Rule induction based on RST ####
 ## determine feature subset/reduct 	
 reduct <- FS.permutation.heuristic.reduct.RST(decision.table,  permutation = NULL)

 ## generate rules						 
 rules.rst <- RI.indiscernibilityBasedRules.RST(decision.table, reduct)
 summary(rules.rst)
 
 ## predicting newdata
 ## in this case, we are using the same dataset as training data
 res.1 <- predict(rules.rst, decision.table[, -ncol(decision.table)])
 print(res.1)