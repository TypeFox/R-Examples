 library(RoughSets)
 data(RoughSetData)
 decision.table <- RoughSetData$hiring.dt 

 ## evaluate single reduct
 res.1 <- FS.quickreduct.RST(decision.table)
 print(res.1)
 
 ## generate new decision table according to the reduct
 new.decTable <- SF.applyDecTable(decision.table, res.1)
 print(new.decTable)