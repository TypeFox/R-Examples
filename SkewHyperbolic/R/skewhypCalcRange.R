######Calculate Range function################################################
skewhypCalcRange <- function(mu = 0, delta = 1, beta = 1, nu = 1,
                             param = c(mu,delta,beta,nu),
                             density = TRUE, tol = 10^(-5), ...){

    ## check parameters
    parResult <- skewhypCheckPars(param)
    case <- parResult$case
    errMessage <- parResult$errMessage
    if(case == "error") stop(errMessage)
    mu <- param[1]
    delta <- param[2]
    beta <- param[3]
    nu <- param[4]
    distMode <- skewhypMode(param = param)

    if (density == TRUE){
        xHigh <- distMode + delta
        while (dskewhyp(xHigh, param = param) > tol){
            xHigh <- xHigh +
                skewhypStepSize(dskewhyp(xHigh, param = param) - tol,
                                delta, beta, nu, side = "right")
        }

        xLow <- distMode  - delta
        while (dskewhyp(xLow, param = param) > tol){
            xLow <- xLow -
                skewhypStepSize(dskewhyp(xLow, param = param) - tol,
                                delta, beta, nu, side = "left")
        }

        ## find xLower and xUpper
        zeroFun <- function(x) dskewhyp(x, param = param) - tol
        xUpper <- uniroot(zeroFun, c(distMode,xHigh),...)$root
        xLower <- uniroot(zeroFun, c(xLow, distMode),...)$root

    } else { # density == FALSE, use probability
        upperProb <- function(x){
            px <- integrate(dskewhyp, x, Inf, param = param)$value
            return(px)
        }

        ## find xHigh, xLow
        xHigh <- distMode + delta
        while (upperProb(xHigh) > tol){
            xHigh <- xHigh +
                skewhypStepSize(upperProb(xHigh) - tol,
                                delta, beta, nu, side = "right")
        }

        lowerProb <- function(x){
            px <- integrate(dskewhyp, -Inf, x, param = param)$value
            return(px)
        }

        xLow <- distMode  - delta
        while (lowerProb(xLow) > tol){
            xLow <- xLow -
                skewhypStepSize(lowerProb(xLow) - tol,
                                delta, beta, nu, side = "left")
        }

        ## find xLower and xUpper
        zeroFun <- function(x) upperProb(x) - tol
        xUpper <- uniroot(zeroFun, c(distMode,xHigh),...)$root
        zeroFun <- function(x) lowerProb(x) - tol
        xLower <- uniroot(zeroFun, c(xLow, distMode),...)$root
    }

    ## put it all together
    range <- c(xLower, xUpper)

    return(range)

}


skewhypStepSize <- function(dist, delta, beta, nu, side = c("right","left"))
{
    ## Purpose: determine the step size for a skewhyperbolic
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: David Scott, Date: 17 Mar 2010, 21:50
    side <- match.arg(side)
    if (beta > 0){
        step <- ifelse(side == "left", delta,
                       delta*abs(beta)*(nu*dist)^(-2/nu))
    }
    if (beta < 0){
        step <- ifelse(side == "right", delta,
                       delta*abs(beta)*(nu*dist)^(-2/nu))
    }
    if (isTRUE(all.equal(beta, 0))){
        step <- exp(dist/nu)
    }
    return(step)
}


