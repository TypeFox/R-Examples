#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[getFpTransforms.R] by DSB Don 04/09/2008 16:49 (CEST) on daniel@puc.home>
##
## Description:
## Transform a variable according to FP transformation formula and attach proper
## names to the resulting design matrix.
##
## History:
## 04/07/2008   copy from thesis function collection.
## 04/09/2008   added argument center
#####################################################################################


getFpTransforms <- function (vec,       # positive (== already shifted and scaled) column vector with proper colname
                             powers,    # power vector with at least one element
                             center=TRUE # center the columns around zero?
                             )
{
    name <- colnames (vec)
    vec <- as.vector (vec)
    len <- length (vec)
    np <- length (powers)

    ret <- matrix (nrow = len, ncol = np)
    retColnames <- character (np)

    lastCol <- 1                        # begin of recursion
    lastPow <- 0
    numRepPow <- 0                      # invariant: already numRepPow times was this power repeated

    for (i in seq_along (powers)){      # invariant: about to write ith column
        if ((pi <- powers[i]) == lastPow && i != 1){ # repeated powers case
            lastCol <- lastCol * log (vec)
            numRepPow <- numRepPow + 1
            retColnames[i] <-           # name is a bit complicated
                if (pi == 0){           # log was repeated
                    paste (getTransformName (name, 0), numRepPow + 1, sep = "^")
                } else {                # other power was repeated
                    tmp <- paste (getTransformName (name, pi), getTransformName (name, 0), sep =  "*")
                    if (numRepPow > 1)
                        tmp <- paste (tmp, numRepPow, sep = "^")
                    tmp
                }
        } else {                        # normal case
            lastCol <- vec %bt% pi
            retColnames[i] <- getTransformName (name, pi)
            numRepPow <- 0
            lastPow <- pi
        }
        ret[, i] <- lastCol
    }

    if(center)
        ret <- scale(ret, center=TRUE, scale=FALSE)
    
    colnames (ret) <- retColnames
    return (ret)
}

####################################################################################################

'%bt%' <- function (x,            # Box Tidwell transformation of vector x (not vectorized for pow!)
                    pow           # by power pow
                    )
{
    if (pow){
        x^pow
    } else {
        log (x)
    }
}

####################################################################################################

getTransformName <- function (name, pow)# helper function
{
    if (pow){
        paste (name, pow, sep = "^")
    } else {
        paste ("log(", name, ")", sep = "")
    }
}
