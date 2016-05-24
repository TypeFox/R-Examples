crimCV <-
function (Dat, ng, dpolyp = 3, dpolyl = 3, model = "ZIPt", rcv = FALSE, 
    init = 20, Risk = NULL) 
{
    no <- ncol(Dat)
    oyr <- 0:(no - 1)/(no - 1)
    if (model == "ZIPt") {
        X <- bs(oyr, degree = dpolyp, intercept = TRUE)
        out <- dmZIPt(Dat, X, ng, rcv, init, Risk)
        return(out)
    }
    else if (model == "ZIP") {
        X <- bs(oyr, degree = dpolyp, intercept = TRUE)
        Z <- bs(oyr, degree = dpolyl, intercept = TRUE)
        out <- dmZIP(Dat, X, Z, ng, rcv, init, Risk)
        return(out)
    }
    else {
        stop("crimCV: The model specified is not currently supported!")
    }
}
