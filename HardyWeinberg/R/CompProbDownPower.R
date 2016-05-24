CompProbDownPower <-
function (nAA, nBB, EnAB, prob, theta, vec = NULL) 
{
    pr <- prob * EnAB * (EnAB - 1)/(theta * (nAA + 1) * (nBB + 1))
    nvec <- c(vec, pr)
    if (EnAB > 3) {
        nvec <- CompProbDownPower(nAA + 1, nBB + 1, EnAB - 2, pr, theta, 
            nvec)
    }
    return(nvec)
}
