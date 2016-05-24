oneSdRule<- function (CVout) 
{
#Changjiang Xu. October 29, 2009
#CV[,1] - cv estimates
#CV[,2] - sd of estimates
    cverrs <- CVout[, 1]
    indMin <- which.min(cverrs)
    fmin <- CVout[indMin,2]
    cutOff <- fmin + cverrs[indMin]
    min(which(cverrs<cutOff)) 
}
