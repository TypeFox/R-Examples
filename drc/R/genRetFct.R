"genRetFct" <- function(fct, parmVec, notFixed)
{
#    ## Defining the model function adjusted for scaling
#    retFct <- function(doseScaling, respScaling, lenData)
#    {   
#        parmMat <- matrix(parmVec / c(1, respScaling, respScaling, doseScaling, 1), lenData, numParm, byrow = TRUE)
#        
#        fct <- function(dose, parm) 
#        {        
#            parmMat[, notFixed] <- parm        
#            cParm <- parmMat[, 2]
#            cParm + (parmMat[, 3] - cParm)/((1+exp(parmMat[, 1]*(log(dose/parmMat[, 4]))))^parmMat[, 5])
#        }
#        fct        
#    }
#    retFct
}