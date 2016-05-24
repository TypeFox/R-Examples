## Helper functions for mean functions
"pickParm" <- function(parmVec, indexVec, parmNo)
{
    function(parm)
    {
        parmVec[indexVec] <- parm
        parmVec[parmNo]
    }
}

"monoParm" <- function(parmVec, indexVec, parmNo, signVal)
{
    function(parm)
    {
        parmVec[indexVec] <- parm
        signVal * parmVec[parmNo]
    }
}
