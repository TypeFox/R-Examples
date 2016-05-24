  library(RoughSets)
  
  dt.ex1 <- data.frame(c(1,0,2,1,1,2,2,0), c(0, 1,0, 1,0,2,1,1), 
                         c(2,1,0,0,2,0,1,1), c(2,1,1,2,0,1,1,0), c(0,2,1,2,1,1,2,1))
 colnames(dt.ex1) <- c("aa", "bb", "cc", "dd", "ee")
 decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 5, indx.nominal = c(1:5))

 P <- c(2,3)
 
 ####### Perform indiscernibility relation #######
 IND <- BC.IND.relation.RST(decision.table, feature.set = P)

 ####### Perform lower and upper approximatino #####
 roughset <- BC.LU.approximation.RST(decision.table, IND)
 
 ####### Determine the positive region ######
 region <- BC.positive.reg.RST(decision.table, roughset) 
 
 print("Indiscernibility Relation:")
 print(IND)
 
 print("Lower and Upper Approximations:")
 print(roughset)
 
 print("Regions:")
 print(region)