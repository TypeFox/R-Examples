##############################################
## Example: Missing Value Completion
##############################################
## Deletion cases
 dt.ex1 <- data.frame(
      c(100.2, 102.6, NA, 99.6, 99.8, 96.4, 96.6, NA), 
      c(NA, "yes", "no", "yes", NA, "yes", "no", "yes"), 
      c("no", "yes", "no", "yes", "yes", "no", "yes", NA),
      c("yes", "yes", "no", "yes", "no", "no", "no", "yes"))
 colnames(dt.ex1) <- c("Temp", "Headache", "Nausea", "Flu")
 decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 4, 
                                     indx.nominal = c(1:4))
 indx1 = MV.deletionCases(decision.table)
 new.decTable1 <- SF.applyDecTable(decision.table, indx1)
 
 ## The most common value
 indx2 = MV.mostCommonValResConcept(decision.table)
 new.decTable2 <- SF.applyDecTable(decision.table, indx2)

 ## Replacing missing attribute values
 ## by the attribute mean/common values
 indx3 = MV.mostCommonVal(decision.table)
 new.decTable3 <- SF.applyDecTable(decision.table, indx3)

 ## Global Closest Fit
 indx4 = MV.globalClosestFit(decision.table)
 new.decTable4 <- SF.applyDecTable(decision.table, indx4)

 ## Concept Closest Fit
 indx5 = MV.conceptClosestFit(decision.table)
 new.decTable5 <- SF.applyDecTable(decision.table, indx5)