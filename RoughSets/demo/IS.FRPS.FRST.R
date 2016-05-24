 #############################################
 ## Example: Evaluate instances/objects and
 ## generate new decision table
 #############################################
 dt.ex1 <- data.frame(c(0.5, 0.2, 0.3, 0.7, 0.2, 0.2), c(0.1, 0.4, 0.2, 0.8, 0.4, 0.4), 
                       c(0, 0, 0, 1, 1, 1))
 colnames(dt.ex1) <- c("a1", "a2", "d")
 decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 3, indx.nominal = c(3))

 ## evaluate instances
 res.1 <- IS.FRPS.FRST(decision.table, type.alpha = "FRPS.3")
 print(res.1)
 
 ## generate new decision table
 new.decTable <- SF.applyDecTable(decision.table, res.1)
 print(new.decTable)