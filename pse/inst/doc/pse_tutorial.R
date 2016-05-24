### R code from vignette source 'pse_tutorial.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: seed
###################################################
set.seed(17)


###################################################
### code chunk number 2: params
###################################################
factors <- c("r", "K", "X0")
q <- c("qnorm", "qnorm", "qunif")
q.arg <- list( list(mean=1.7, sd=0.3), list(mean=40, sd=1), 
	list(min=1, max=50) )


###################################################
### code chunk number 3: qdata
###################################################
qdata <- function(p, data) quantile(x=data, probs=p)


###################################################
### code chunk number 4: qdunif
###################################################
qdunif<-function(p, min, max) floor(qunif(p, min, max))


###################################################
### code chunk number 5: model
###################################################
oneRun <- function (r, K, Xo) {
    X <- Xo
    for (i in 0:20) {
       X <- X+r*X*(1-X/K)
    }   
    return (X) 
}
modelRun <- function (my.data) {
	return(mapply(oneRun, my.data[,1], my.data[,2], my.data[,3]))
}


###################################################
### code chunk number 6: LHS
###################################################
library(pse)
myLHS <- LHS(modelRun, factors, 200, q, q.arg, nboot=50)


###################################################
### code chunk number 7: ecdf
###################################################
plotecdf(myLHS)


###################################################
### code chunk number 8: corplot
###################################################
plotscatter(myLHS)


###################################################
### code chunk number 9: prcc
###################################################
plotprcc(myLHS)


###################################################
### code chunk number 10: pic
###################################################
pic(myLHS, nboot=40)


###################################################
### code chunk number 11: sbma
###################################################
newLHS <- LHS(modelRun, factors, 250, q, q.arg)
(mySbma <- sbma(myLHS, newLHS))


###################################################
### code chunk number 12: p2
###################################################
factors <- c("r", "K", "X0")
q <- c("qnorm", "qnorm", "qunif")
q.arg <- list( list(mean=1.7, sd=0.3), list(mean=40, sd=1), 
	list(min=1, max=50) )
Time <- 6
oneRun <- function (r, K, Xo) {
	X <- array();
	X[1] <- Xo; # Caution, X1 gets overwritten
	for (i in 1:Time) {
		Xl <- X[length(X)]
		X[i] <- Xl + r*Xl*(1-Xl/K)
	}
	return (X)
}
modelRun <- function (dados) {
	mapply(oneRun, dados[,1], dados[,2], dados[,3])
}


###################################################
### code chunk number 13: multiLHS
###################################################
res.names <- paste("Time",1:Time)
myLHS <- LHS(modelRun, factors, 100, q, q.arg, res.names, nboot=50)


###################################################
### code chunk number 14: ecdf2
###################################################
plotecdf(myLHS, stack=TRUE)


###################################################
### code chunk number 15: corplot2
###################################################
plotscatter(myLHS, index.res=c(1,3,6), add.lm=FALSE)


###################################################
### code chunk number 16: prcc2
###################################################
plotprcc(myLHS, index.res=c(1,3,6))


###################################################
### code chunk number 17: target
###################################################
targetLHS <- target.sbma (target=0.3, modelRun, factors, 
	q, q.arg, res.names, FUN=min) 


###################################################
### code chunk number 18: pse_tutorial.Rnw:432-434
###################################################
uncoupledLHS <- LHS(model=NULL, factors, 50, q, q.arg)
write.csv(get.data(uncoupledLHS), file="mydata.csv")


###################################################
### code chunk number 19: pse_tutorial.Rnw:443-444
###################################################
myresults <- apply(get.data(uncoupledLHS), 1, mean)


###################################################
### code chunk number 20: pse_tutorial.Rnw:447-448
###################################################
coupledLHS <- tell(uncoupledLHS, myresults)


