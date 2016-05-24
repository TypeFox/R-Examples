 library(RoughSets)
 
 ###########################################################
 ##### 1. Example: Using simple decision table containing 
 #####             nominal values on decision attribute
 ###########################################################
 dt.ex1 <- data.frame(c(-0.4, -0.4, -0.3, 0.3, 0.2, 0.2), 
                      c(-0.3, 0.2, -0.4, -0.3, -0.3, 0),
				        c(-0.5, -0.1, -0.3, 0, 0, 0),
				        c("no", "yes", "no", "yes", "yes", "no"))
 colnames(dt.ex1) <- c("a", "b", "c", "d")
 decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 4)

 ## let us consider the first and second attributes only as conditional attribute
 condAttr <- c(1, 2)
 
 ## let us consider the fourth attribute as decision attribute
 decAttr <- c(4)
 
 #### Calculate fuzzy indiscernibility relation ####
 control.ind <- list(type.aggregation = c("t.tnorm", "lukasiewicz"), 
                     type.relation = c("tolerance", "eq.1"))
 control.dec <- list(type.aggregation = c("crisp"), type.relation = "crisp")

 IND.condAttr <- BC.IND.relation.FRST(decision.table, attributes = condAttr, 
                                      control = control.ind) 
 IND.decAttr <- BC.IND.relation.FRST(decision.table, attributes = decAttr, 
                                      control = control.dec) 
 
 #### Calculate fuzzy lower and upper approximation using type.LU : 
 #### "implicator.tnorm" 
 control <- list(t.implicator = "lukasiewicz")
 FRST.LU <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
               type.LU = "implicator.tnorm", control = control)

 #### Determine fuzzy regions ####
 fuzzy.region <- BC.positive.reg.FRST(decision.table, FRST.LU)
 
 
 print("Indiscernibility Relation: ")
 print(IND.condAttr)
 print(IND.decAttr)
 
 print("Fuzzy Lower and Upper Approximations:")
 print(FRST.LU)
 
 print("Fuzzy Region:")
 print(fuzzy.region)