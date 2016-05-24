predict.dpa1 <-
function (object, newdata, ...) 
{
    if(is.matrix(newdata))
        nbreDExempleMax=dim(newdata)[1]
    else
        nbreDExempleMax=1
    res=c()
    for(nbreDExemple in 1:nbreDExempleMax)
    {
        maximum = -1
        maximumIndice = -1
        supposonsLaclCle = 0
        for (dataCleI in object$mean) {
            if(is.matrix(newdata))
               corr <- cor(newdata[nbreDExemple,], dataCleI)
            else
               corr <- cor(newdata, dataCleI)
            if (maximum < corr) {
               maximumIndice = supposonsLaclCle
               maximum = corr
            }
            supposonsLaclCle= supposonsLaclCle+1
        }
        res = c(res,maximumIndice)
    }
    return (factor(res))
}