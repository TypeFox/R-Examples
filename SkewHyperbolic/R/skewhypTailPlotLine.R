#----------------------------------------------------------

skewhypTailPlotLine <- function(x, mu = 0, delta = 1, beta = 1, nu = 1,
                                param = c(mu,delta,beta,nu),
                                side = c("right", "left"), ...)
{
    ## Purpose: Add skew hyperbolic distribution line to tail plot
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: David Scott, Date: 31 Mar 2010, 23:38

    ## check parameters
    parResult <- skewhypCheckPars(param)
    case <- parResult$case
    errMessage <- parResult$errMessage
    if (case == "error")
        stop(errMessage)
    param <- as.numeric(param)
    mu <- param[1]
    delta <- param[2]
    beta <- param[3]
    nu <- param[4]
    side <- match.arg(side)

    x <- sort(as.numeric(x))
    if (side =="left"){
        lines(x, pskewhyp(x, param = param), ...)
    }
    if (side =="right"){
        lines(x, 1 - pskewhyp(x, param = param), ...)
    }
    invisible(NULL)
}

