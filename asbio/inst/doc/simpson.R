### R code from vignette source 'simpson.rnw'

###################################################
### code chunk number 1: closeConnetions
###################################################
allCon <- showConnections()
socketCon <- as.integer(rownames(allCon)[allCon[, "class"] == "sockconn"])
sapply(socketCon, function(ii) close.connection(getConnection(ii)) )


