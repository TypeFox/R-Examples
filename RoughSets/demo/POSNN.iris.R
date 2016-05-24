library(RoughSets)
 #############################################################
 ## In this example, we are using Iris dataset.
 ## It should be noted that since values of decision attribute in string,
 ## so they should be transformed into numeric values using unclass()
 #############################################################
 data(iris)
 set.seed(2)
 
 irisShuffled <- iris[sample(nrow(iris)),]
 irisShuffled[,5] <- unclass(irisShuffled[,5])
 iris.training <- irisShuffled[1:105,]
 real.iris <- matrix(irisShuffled[106:nrow(irisShuffled),5], ncol = 1)
 colnames(iris.training) <- c("Sepal.Length", "Sepal.Width", "Petal.Length", 
                        "Petal.Width", "Species")
 decision.table <- SF.asDecisionTable(dataset = iris.training, decision.attr = 5, indx.nominal = c(5))
		
 ## define newdata
 tst.iris <- SF.asDecisionTable(dataset = irisShuffled[106:nrow(irisShuffled),1:4])
	   
 ## perform FRNN algorithm using lower/upper approximation: Implicator/tnorm based approach
 control <- list(type.LU = "implicator.tnorm", k = 20, t.tnorm = "lukasiewicz", 
                 type.relation = c("tolerance", "eq.1"), t.implicator = "lukasiewicz")

 res.test.POSNN <- C.POSNN.FRST(decision.table = decision.table, 
                            newdata = tst.iris, control = control)

 err = 100*sum(real.iris!=res.test.POSNN)/nrow(real.iris)

print("The result: ")
print(res.test.POSNN)
print("POSNN: percentage Error on Iris: ")
print(err) 						   							