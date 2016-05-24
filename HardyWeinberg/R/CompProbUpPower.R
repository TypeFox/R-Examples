CompProbUpPower <-
function (nAA, nBB, EnAB, prob, MaxHet, theta, vec = NULL) 
{
    pr <- prob * theta * nAA * nBB/((EnAB + 2) * (EnAB + 1))
    nvec <- c(vec, pr)
    if (EnAB < MaxHet - 2) {
        nvec <- CompProbUpPower(nAA - 1, nBB - 1, EnAB + 2, pr, MaxHet, theta,  
            nvec)
    }
    return(nvec)
}
