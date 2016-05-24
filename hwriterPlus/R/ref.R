### Convenience functions to retrieve a reference number, given a label
###
### DJS, 19/1/2012

eqRef <- function(label){
    eqLabel <- paste("eq:", label, sep = "")
    eqnList <- get("hwriterEquationList", .hwriterGlobalEnv)
    eqNum <- which(eqnList == eqLabel)
    eqNum
}

