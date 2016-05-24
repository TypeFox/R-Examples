mci.transmat <-
function (mcidataset, submarkets, suppliers, mcivariable1, ...) {
mcidataset_rows <- nrow(mcidataset)
addmcivars <- unlist(list(...))
addmcivars_count <- length(addmcivars)
mcivariablelog1 <- mci.transvar (mcidataset, submarkets, suppliers, mcivariable1, output_ij=TRUE)
v <- 0
addmcivariablelog <- data.frame(matrix(0, nrow=mcidataset_rows, ncol=addmcivars_count))
varname <- character()

for (v in 1:addmcivars_count) {
varname <- addmcivars[[v]]
addmcivariablelog[v] <- mci.transvar (mcidataset, submarkets, suppliers, varname, output_ij=FALSE)
colnames(addmcivariablelog)[v] <- paste(varname, "_t", sep="")
}
mcilinoutput <- data.frame(mcivariablelog1, addmcivariablelog)
return(mcilinoutput)
}
