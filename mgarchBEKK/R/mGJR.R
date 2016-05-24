##' Bivariate GJR Estimation
##'
##' Provides bivariate GJR (\code{mGJR(p,q,g)}) estimation procedure.
##'
##'
##'
##' @param eps1 First time series.
##' @param eps2 Second time series.
##' @param order mGJR(p, q, g) order a three element integer vector
##'     giving the order of the model to be fitted. \code{order[2]}
##'     refers to the ARCH order and \code{order[1]} to the GARCH
##'     order and \code{order[3]} to the GJR order.
##' @param params Initial parameters for the \code{optim} function.
##' @param fixed A two dimensional vector that contains the user
##'     specified fixed parameter values.
##' @param method The method that will be used by the \code{optim}
##'     function. See \code{?optim} for available options.
##' @return Estimation results packaged as \code{mGJR} class instance. The values are defined as:
##' \describe{
##'    \item{eps1}{first time series}
##'    \item{eps2}{second time series}
##'    \item{length}{length of each series}
##'    \item{order}{order of the mGJR model fitted}
##'    \item{estimation.time}{time to complete the estimation process}
##'    \item{total.time}{time to complete the whole routine within the mGJR.est process}
##'    \item{estimation}{estimation object returned from the optimization process, using \code{optim}}
##'    \item{aic}{the AIC value of the fitted model}
##'    \item{est.params}{estimated parameter matrices}
##'    \item{asy.se.coef}{asymptotic theory estimates of standard errors of estimated parameters}
##'    \item{cor}{estimated conditional correlation series}
##'    \item{sd1}{first estimated conditional standard deviation series}
##'    \item{sd2}{second estimated conditional standard deviation series}
##'    \item{H.estimated}{estimated series of covariance matrices}
##'    \item{eigenvalues}{estimated eigenvalues for sum of Kronecker products}
##'    \item{uncond.cov.matrix}{estimated unconditional covariance matrix}
##'    \item{resid1}{first estimated series of residuals}
##'    \item{resid2}{second estimated series of residuals}
##' }
##'
##' @references{
##'   Bauwens L., S. Laurent, J.V.K. Rombouts, Multivariate GARCH models: A survey, April, 2003
##'
##'   Bollerslev T., Modelling the coherence in short-run nominal exchange rate: A multivariate generalized ARCH approach, Review of Economics and Statistics, 498--505, 72, 1990
##'
##'   Engle R.F., K.F. Kroner, Multivariate simultaneous generalized ARCH, Econometric Theory, 122-150, 1995
##'
##'   Engle R.F., Dynamic conditional correlation: A new simple class of multivariate GARCH models, Journal of Business and Economic Statistics, 339--350, 20, 2002
##'
##'   Tse Y.K., A.K.C. Tsui, A multivariate generalized autoregressive conditional heteroscedasticity model with time-varying correlations, Journal of Business and Economic Statistics, 351-362, 20, 2002
##' }
##'
##'
##' @examples
##' \dontrun{
##'   sim = BEKK.sim(1000)
##'   est = mGJR(sim$eps1, sim$eps2)
##' }
##'
##' @import stats
##' @useDynLib mgarchBEKK
##' @export
mGJR <- function(eps1, eps2, order=c(1,1,1), params=NULL, fixed=NULL, method="BFGS") {
    ## check the given time series
    if(length(eps1) != length(eps2))
    {
        stop("time series are different in length")
    }

    ## check the given order
    ## orders should be integers
    if(order[1] != as.integer(order[1]) || order[2] != as.integer(order[2]) || order[2] != as.integer(order[2]))
    {
        stop("order should contain integer values")
    }

    ## GARCH and GJR effect could be set to 0, but, ARCH could not be 0
    if(order[1] < 0 || order[3] < 0 || order[2] < 1)
    {
        stop("mGJR(",order[1],",",order[2],",",order[3],") is not implemented.")
    }

    ## construct the paramters list.
    if(order[3] > 0)
    {
        tempw = 1
    }
    else
    {
        tempw = 0
    }

    length.params = 3 + (order[1] * 4) + (order[2] * 4) + (order[3] * 4) + tempw # set the length of the parameter list

    if(is.null(params))
    {
        ## WARNING
        ## for order = 0x1x0 1x1x0, 1x2x0, 2x1x0,2x2x0 we do offer some initial parameter lists.
        ## for other trials, the remaining parameters are set to `0'
        if(order[1] == 0 && order[2] == 1 && order[3] == 0)
        {
            params = c(2, 0, 2, 0.4, 0.1, 0.1, 0.4)
        }
        else if(order[1] == 1 && order[2] == 1 && order[3] == 0)
        {
            params = c(2, 0, 2, 0.4, 0.1, 0.1, 0.4, 0.4, 0.1, 0.1, 0.4)
        }
        else if(order[1] == 2 && order[2] == 1 && order[3] == 0)
        {
            params = c(2, 0, 2, 0.4, 0.1, 0.1, 0.4, 0.4, 0.1, 0.1, 0.4, 0.2, 0.1, 0.1, 0.2)
        }
        else if(order[1] == 1 && order[2] == 2 && order[3] == 0)
        {
            params = c(2, 0, 2, 0.4, 0.1, 0.1, 0.4, 0.2, 0.1, 0.1, 0.2, 0.4, 0.1, 0.1, 0.4)
        }
        else if(order[1] >= 2 && order[2] >= 2 && order[3] == 0)
        {
            params = c(2, 0, 2, 0.4, 0.1, 0.1, 0.4, 0.2, 0.1, 0.1, 0.2, rep(0,(order[2] - 2) * 4), 0.4, 0.1, 0.1, 0.4, 0.2, 0.1, 0.1, 0.2, rep(0,(order[1] - 2) * 4))
        }
        else if(order[1] >= 1 && order[2] >= 1 && order[3] >= 0)
        {
            params = c(2, 0, 2, 0.4, 0.1, 0.1, 0.4, rep(0,(order[2] - 1) * 4), 0.4, 0.1, 0.1, 0.4, rep(0,(order[1] - 1) * 4), rep(0.1,(order[3]) * 4), 0.5)
        }
        else
        {
            params = c(rep(0, length.params - 1), 0.5)
        }

        cat("\nWarning: initial values for the parameters are set at:\n\t", params,"\n")
    }
    else if(length(params) != length.params)
    {
        stop("length of the initial parameter list doesn't conform required length (3 + (order[1] * 4) + (order[2] * 4) + (order[3] * 4) + 1).");
    }

    ## check the given fixed parameters
    if(!is.null(fixed))
    {
        ## check the format of the fixed parameters
        if(!is.array(fixed))
        {
            stop("fixed should be an array of two vectors. Try fixed = array(c(a,b,c,d,...), dim = c(2,y))")
        }

        if(dim(fixed)[1] != 2)
        {
            stop("fixed should be an array of two vectors. Try fixed = array(c(a,b,c,d,...), dim = c(2,y))")
        }

        if(length(fixed[1,]) != length(fixed[2,]))
        {
            stop("fixed should be an array of two vectors. Try fixed = array(c(a,b,c,d,...), dim = c(2,y))")
        }

        ## check the first dimension, if it contains appropriate index values,
        ## that is integer values rather than floating or negative numbers
        for(count in 1:length(fixed[1,]))
        {
            if((fixed[1,count] != as.integer(fixed[1,count])) || (fixed[1,count] <= 0))
            {
                stop("First dimension of the fixed array should contain only positive integer values for indexing purposes")
            }
        }

        ## check the length of the fixed parameters
        if(length(fixed[1,]) > length(params))
        {
            stop("fixed array could not contain more index-value pairs than the params array length");
        }
    }

    ## check the method specified in the argument list
    if(!(
        (method == "Nelder-Mead") ||
        (method == "BFGS") ||
        (method == "CG") ||
        (method == "L-BFGS-B") ||
        (method == "SANN")
    ))
    {
        stop("'", method, "' method is not available")
    }


    fake.params = params
    if(!is.null(fixed))
    {
        ## extract the parameters specified in the fixed list.
        fake.params = params
        for(i in 1:length(fixed))
        {
            fake.params[fixed[1,][i]] = NA
        }
        fake.params = na.omit(fake.params)
    }

    ## parameters seem appropriate
    ## define the loglikelihood function
    loglikelihood.GJR <- function(params)
    {
        loglikelihood.GJR <- .C("loglikelihood_GJR",
                                as.vector(params, mode = "double"),
                                as.vector(fixed[1,], mode = "integer"),
                                as.vector(fixed[2,], mode = "double"),
                                as.integer(length(fixed[1,])),
                                as.vector(eps1, mode = "double"),
                                as.vector(eps2, mode = "double"),
                                as.integer(length(eps1)),
                                as.vector(order, mode = "integer"),
                                retval = 0.0,
                                PACKAGE = "mgarchBEKK"
                                )

        if(is.nan(loglikelihood.GJR$retval) == T)
        {
            nonusedret = 1e+100

        }
        else
        {
            nonusedret = loglikelihood.GJR$retval
        }
        nonusedret
    }

    ## begin estimation process

    ## first log the start time
    start = Sys.time()
    cat("Starting estimation process via loglikelihood function implemented in C.\n")
    cat("Optimization Method is '", method, "'\n")

    ## call the optim function
    estimation = optim(fake.params, loglikelihood.GJR, method = method, hessian = T)
    ## estimation completed

    cat("Estimation process completed.\n")


    ## log estimation time
    est.time = difftime(Sys.time(), start)

    ## calculate the AIC
    ## it is estimation value + number of estimated parameters
    aic = estimation$value + (length(params) - length(fixed[1,]))

    ## following script will prepare an object that holds the estimated
    ## parameters and some useful diagnostics data like estimated correlation,
    ## standard deviation, eigenvalues and so on.

    ## TODO
    ## estimation$hessian is non-existing if fixed parameter list contains all the
    ## paramters to be estimated. That is that the estimation procedure gets no parameters,
    ## thus, there is no errors... Fix it... How? Whether encapsulate with an "if" statement, probably not
    ## efficient, or give a fake hessian
    ## give a fake hessian
    if(length(fake.params) == 0)
    {
        estimation$hessian = matrix(c(0,0.1,0.2,0), nrow = 2, ncol = 2)
    }

                                        # get the hessian matrix and grap the diagonal
    inv.hessian.mat = solve(estimation$hessian)

    diag.inv.hessian = array(rep(1,dim(inv.hessian.mat)[1]))
    for(count in 1:dim(inv.hessian.mat)[1])
    {
        diag.inv.hessian[count] = sqrt(inv.hessian.mat[count,count])
    }

    ## fix the asymptotic-theory standard errors of the
    ## coefficient estimates with fixed parameters
    if(!is.null(fixed))
    {
        temp.diag.inv.hessian = numeric()
        shifted = 0
        for(count in 1:length.params)
        {
            check.point = 0
            for(i in 1:length(fixed[1,]))
            {
                if(count == fixed[1,i])
                {
                    check.point = 1
                    shifted = shifted + 1
                    temp.diag.inv.hessian[count] = 0
                    break
                }
            }
            if(check.point == 0)
            {
                temp.diag.inv.hessian[count] = diag.inv.hessian[count - shifted]
            }
        }
        diag.inv.hessian = temp.diag.inv.hessian
    }

    ## construct the asymptotic-theory standard errors of the coefficient estimates matrices
    parnum = 1 + order[1] + order[2] + order[3]    # calculate number of paramater matrices
    asy.se.coef = list()                # declare the asy.se.coef matrices list

    ## first initialize the first asy.se.coef matrix, corresponding to the C matrix
    asy.se.coef[[1]] = array(c(diag.inv.hessian[1], 0, diag.inv.hessian[2], diag.inv.hessian[3]), dim = c(2,2))

    ## following loop initalizes the ARCH and GARCH parameter matrices respectively
    for(count in 1:(parnum - 1))
    {
        asy.se.coef[[count + 1]] = array(diag.inv.hessian[(4 + (count - 1) * 4):(8 + (count - 1) * 4)], dim = c(2,2));
    }

    asy.se.coef[[parnum + 1]] = diag.inv.hessian[length.params]

    buff.par = list()                    # declare the parameter list

    ## shift the fixed parameters inside the estimated paramters
    if(!is.null(fixed))
    {
        estim.params = numeric()
        shifted = 0
        for(count in 1:length.params)
        {
            check.point = 0
            for(i in 1:length(fixed[1,]))
            {
                if(count == fixed[1,i])
                {
                    check.point = 1
                    shifted = shifted + 1
                    estim.params[count] = fixed[2,i]
                    break
                }
            }
            if(check.point == 0)
            {
                estim.params[count] = estimation$par[count - shifted]
            }
        }
    }
    else
    {
        estim.params = estimation$par
    }

    ## first initialize the C matrix
    buff.par[[1]] = array(c(estim.params[1], 0, estim.params[2], estim.params[3]), dim = c(2,2))

    ## following loop initalizes the ARCH and GARCH parameter matrices respectively
    for(count in 1:(parnum - 1))
    {
        buff.par[[count + 1]] = array(estim.params[(4 + (count - 1) * 4):(8 + (count - 1) * 4)], dim = c(2,2));
    }

    buff.par[[parnum + 1]] = estim.params[length.params]

    ## calculate the transposes of the parameter matrices
    buff.par.transposed = list()
    for(count in 1:parnum)
    {
        buff.par.transposed[[count]] = t(buff.par[[count]])
    }

    ## start diagnostics
    cat("Starting diagnostics...\n")
    cat("Calculating estimated:\n 1. residuals,\n 2. correlations,\n 3. standard deviations,\n 4. eigenvalues.\n")

    HLAGS = list()        # list of H lags that will be used later in the MGARCH implementation
    for(count in 1:max(order))
    {
        HLAGS[[count]] = array(c(1,0,0,1), dim = c(2,2))
    }

    T = length(eps1)    # length of the series

    resid1 = numeric()    # declare the first  residual series
    resid2 = numeric()    # declare the second residual series

    ## initialize the first residuals we are not able to calculate
    for(count in 1:max(order))
    {
        resid1[count] = 0
        resid2[count] = 0
    }

    resid = array(c(0,0), dim = c(2,1)) # declare a temporary residuals buffer

    ## calculate eigenvalues
    temp = 0
    for(count in 2:parnum)
    {
        temp = temp + kronecker(buff.par[[count]], buff.par[[count]])
    }
    eigenvalues = svd(temp)$d

    ## compute the unconditional covariance matrix
    numerat = t(buff.par[[1]]) %*% buff.par[[1]]
    dim(numerat) = c(4,1)
    denom = solve(diag(c(1,1,1,1)) - temp)
    sigma = denom %*% numerat
    dim(sigma) = c(2,2)

    H = sigma    # to initialize: assign the unconditional covariance matrix to the H matrix :)
    H.estimated = array(c(var(eps1), cov(eps1, eps2), cov(eps1, eps2), var(eps2)), dim = c(2,2,T))    # declare the estimated H series
    cor = numeric()        # declare the estimated correlation series
    sd1 = numeric()        # declare the first  estimated standard deviation series
    sd2 = numeric()        # declare the second estimated standard deviation series
    eps = array(c(0,0), dim = c(2,1))    # declare a temporary eps buffer
    CTERM = buff.par.transposed[[1]] %*% buff.par[[1]] # calculate the C'C term

    for(count in (max(order) + 1):T)    # cruical loop! initializing diagnostics data
    {
        ## do the swap calculation for H terms
        if(order[1] >= 2)
        {
            for(tmp.count in max(order):2)
            {
                HLAGS[[tmp.count]] = HLAGS[[(tmp.count - 1)]]
            }
        }
        HLAGS[[1]] = H

        ## a bit complicated but following explanation will be useful hopefully
        ## H = (C')x(C) + (A')(E_t-1)(E_t-1')(A) + (B')(E_t-2)(E_t-2')(B) + ... +  (G')(H_t-1)(G) + (L')(H_t-2)(L) + ...
        ##                    |_____________|          |_____________|             |____________|   |____________| |_____|
        ##                        E1 TERM                  E2 TERM                     G1 TERM         G2 TERM     G3.G4..
        ##                |____________________|   |____________________| |_____|
        ##                        A1 TERM                  A2 TERM        A3.A4..
        ##     |______|  |_____________________________________________________|  |______________________________________|
        ##      C TERM                         A TERM                                              G TERM

        H = CTERM
        ord1 = 1
        for(tmp.count in 1:(order[2] + order[1]))
        {
            if(tmp.count <= order[2])
            {
                ## ARCH EFFECT (A TERM)
                eps = array(c(eps1[count - tmp.count], eps2[count - tmp.count]), dim = c(2,1)) # E TERM
                H = H + buff.par.transposed[[(tmp.count + 1)]] %*% eps %*% t(eps) %*% buff.par[[(tmp.count + 1)]]
            }
            else
            {
                ## GARCH EFFECT (G TERM)
                H = H + buff.par.transposed[[(tmp.count + 1)]] %*% HLAGS[[ord1]]  %*% buff.par[[(tmp.count + 1)]]
                ord1 = ord1 + 1
            }
        }

        for(tmp.count in 1:order[3])
        {
            ## ARCH EFFECT (A TERM)
            eps = array(c(eps1[count - tmp.count], eps2[count - tmp.count]), dim = c(2,1)) # E TERM
            ## THE FOLLOWING MODIFIED BY HARALD, 2004-06-18
            ## if(buff.par[[parnum + 1]]*eps1[count - tmp.count] + (1 - buff.par[[parnum + 1]])*eps2[count - tmp.count] > 0)
            ## {
            ##     S = 1
            ## }
            ## else
            ## {
            ##     S = 0
            ## }
            S = 1-0.5*((cos(pi/4+buff.par[[parnum + 1]])*eps1[count - tmp.count]+sin(pi/4+buff.par[[parnum + 1]])*eps2[count - tmp.count])/sqrt(eps1[count - tmp.count]^2+eps2[count - tmp.count]^2)+1)
            H = H + S*buff.par.transposed[[(tmp.count + 1 + order[1] + order[2])]] %*% eps %*% t(eps) %*% buff.par[[(tmp.count + 1 + order[1] + order[2])]]
        }

        ## TODO add appropriate comments for following assignments and calculations
        H.estimated[,,count] = H
        svdH = svd(H)
        sqrtH = svdH$u %*% diag(sqrt(svdH$d)) %*% t(svdH$v)
        eps = array(c(eps1[count], eps2[count]), dim = c(2,1))

        invsqrtH = solve(sqrtH)
        resid = invsqrtH %*% eps
        resid1[count] = resid[1,1]
        resid2[count] = resid[2,1]

        cor[count] = H[1,2]/(sqrt(H[1,1] * H[2,2]))
        sd1[count] = sqrt(H[1,1])
        sd2[count] = sqrt(H[2,2])
    }


    ## diagnostics ready
    cat("Diagnostics ended...\n")

    names(order) <- c("GARCH component", "ARCH component", "HJR component")
    names(buff.par) <- as.integer(seq(1, parnum + 1))

    retval <- list(
        eps1 = eps1,
        eps2 = eps2,
        series.length = T,
        estimation.time = est.time,
        total.time = difftime(Sys.time(), start),
        order = order,
        estimation = estimation,
        aic = aic,
        asy.se.coef = asy.se.coef,
        est.params = buff.par,
        cor = cor,
        sd1 = sd1,
        sd2 = sd2,
        H.estimated = H.estimated,
        eigenvalues = eigenvalues,
        uncond.cov.matrix = sigma,
        resid1 = resid1,
        resid2 = resid2
    )

    class(retval) = "mGJR"

    cat("Class attributes are ready via following names:\n")
    cat(names(retval), "\n")

    return(retval)
}
