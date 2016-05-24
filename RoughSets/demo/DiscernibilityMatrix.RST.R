 library(RoughSets)
 
 ##########################################################
 ## Example 1: Determining reducts by discernibility matrix
 ##########################################################
 dt.ex1 <- data.frame(c(1,0,2,1,1,2,2,0), c(0, 1,0, 1,0,2,1,1), 
                         c(2,1,0,0,2,0,1,1), c(2,1,1,2,0,1,1,0), c(0,2,1,2,1,1,2,1))
 colnames(dt.ex1) <- c("aa", "bb", "cc", "dd", "ee")
 decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 5, indx.nominal = c(1:5))

 res.1 <- BC.discernibility.mat.RST(decision.table, range.object = NULL)
 