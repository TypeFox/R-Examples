 library(RoughSets)

 data(iris)
 set.seed(2)
 
 irisShuffled <- iris[sample(nrow(iris)),]
 ## transform into numeric values
 irisShuffled[,5] <- unclass(irisShuffled[,5])
 iris.training <- irisShuffled[1:105,]
 real.iris <- matrix(irisShuffled[106:nrow(irisShuffled),5], ncol = 1)
 colnames(iris.training) <- c("Sepal.Length", "Sepal.Width", "Petal.Length", 
                        "Petal.Width", "Species")
 decision.table <- SF.asDecisionTable(dataset = iris.training, decision.attr = 5, indx.nominal = c(5))

 ## Define newdata for testing
 tst.iris <- SF.asDecisionTable(dataset = irisShuffled[106:nrow(irisShuffled),1:4])

 ###### perform FRNN algorithm using lower/upper approximation: 
 ###### Implicator/tnorm based approach
 control <- list(type.LU = "implicator.tnorm", k = 20, 
                 type.aggregation = c("t.tnorm", "lukasiewicz"), 
                 type.relation = c("tolerance", "eq.1"), t.implicator = "lukasiewicz") 									   
 res.1 <- C.FRNN.FRST(decision.table = decision.table, newdata = tst.iris,
                              control = control)

 ###### perform FRNN algorithm using VQRS
 control <- list(type.LU = "vqrs", k = 20, q.some = c(0.1, 0.6), q.most = c(0.2, 1), 
                  type.relation = c("tolerance", "eq.1"), 
                  type.aggregation = c("t.tnorm","lukasiewicz"))
 res.2 <- C.FRNN.FRST(decision.table = decision.table, newdata = tst.iris,
                              control = control)

 ## error calculation
 err.1 = 100*sum(real.iris!=res.1)/nrow(real.iris)
 err.2 = 100*sum(real.iris!=res.2)/nrow(real.iris)

print("FRNN: percentage Error on Iris: ")
print(err.1)
print(err.2) 						   