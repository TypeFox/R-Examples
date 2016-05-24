##' Simulate BEKK processes
##'
##' Provides a procedure to simulate BEKK processes.
##'
##' \code{simulateBEKK} simulates an N dimensional \code{BEKK(p,q)}
##' model for the given length, order list, and initial parameter list
##' where \code{N} is also specified by the user.
##'
##' @param series.count The number of series to be simulated.
##' @param T The length of series to be simulated.
##' @param order BEKK(p, q) order. An integer vector of length 2
##'     giving the orders of the model to fit. \code{order[2]} refers
##'     to the ARCH order and \code{order[1]} to the GARCH order.
##' @param params A vector containing a sequence of parameter
##'     matrices' values.
##' @return Simulated series and auxiliary information packaged as a
##'     \code{simulateBEKK} class instance. Values are:
##'     \describe{
##'         \item{length}{length of the series simulated}
##'         \item{order}{order of the BEKK model}
##'         \item{params}{a vector of the selected parameters}
##'         \item{true.params}{list of parameters in matrix form}
##'         \item{eigenvalues}{computed eigenvalues for sum of Kronecker products}
##'         \item{uncond.cov.matrix}{unconditional covariance matrix of the process}
##'         \item{white.noise}{white noise series used for simulating the process}
##'         \item{eps}{a list of simulated series}
##'         \item{cor}{list of series of conditional correlations}
##'         \item{sd}{list of series of conditional standard deviations}
##'     }
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
##' @examples
##' ## Simulate series:
##' simulated = simulateBEKK(2, 1000, c(1,1))
##'
##' @import stats
##' @export
simulateBEKK <- function(series.count, T, order  = c(1, 1), params = NULL) {
    count.triangular <- function(dimension){
        if(dimension <= 0){
            0
        }
      else{
          dimension + count.triangular(dimension - 1)
      }
    }

    ## check the given order
    ## orders should be integers
    if(order[1] != as.integer(order[1]) || order[2] != as.integer(order[2]))
    {
        stop("order should contain integer values")
    }

    ## GARCH effect could be set to 0, but, ARCH could not be 0
    if(order[1] < 0 || order[2] < 1)
    {
        stop("BEKK(",order[1],",",order[2],") is not implemented.")
    }

    ## init the initial parameters
    params.length = count.triangular(series.count) + (order[1] * series.count^2) + (order[2] * series.count^2)
    if(is.null(params))
    {
        if(order[1] == 1 && order[2] == 1)
        {
            if(series.count == 2)
            {
                params = c(1, 0.8, 1, 0.5, -0.4, 0, 0.3, 0.4, -0.3, 0.5, 0.8)
            }
            else
            {
                params = c(1, 0.2, 1.04, 0.3, 0.01, 0.9,
                           0.3, -0.02, -0.01, 0.01, 0.4, -0.06, 0.02, 0.3, 0.5,
                           0.2, 0.01, -0.1, -0.03, 0.3, -0.06, 0.7, 0.01, 0.5)
            }
        }
        else if(order[1] == 1 && order[2] == 2)
        {
            params = c(1, 0.8, 1, 0.6, -0.4, 0, 0.4, 0.2, 0.1, -0.1, 0.3, 0.6, -0.2, 0.3, 0.5)
        }
        else if(order[1] == 2 && order[2] == 1)
        {
            params = c(1, 0.8, 1, 0.6, -0.4, 0, 0.4, 0.6, -0.2, 0.3, 0.5, 0.2, 0.1, -0.1, 0.3)
        }
        else if(order[1] == 2 && order[2] == 2)
        {
            params = c(1, 0.8, 1, 0.6, -0.4, 0, 0.3, 0.2, 0.1, -0.1, 0.3, 0.6, -0.2, 0.3, 0.5, 0.2, 0.1, -0.1, 0.3)
        }
        else if(order[1] == 0 && order[2] == 1)
        {
            params = c(1, 0.8, 1, 0.9, 0.3, -0.4, 0.8)
        }

        ## fill the rest: might be a bad idea
        params = c(params, rep(0.01, params.length - length(params)))
    }


    ## check the parameter list
    if(length(params) != params.length)
    {
        stop("length of parameters doesn't match the requiered length for the requested BEKK model");
    }

    ## how many parameter matrices in total
    total.par.matrices = 1 + order[1] + order[2]


    ## declare the parameter list
    buff.par = list()

    ## first initialize the C matrix
    tmp.array = array(rep(0, series.count^2), dim = c(series.count, series.count))
    iter = 1
    for(i in 1:series.count)
    {
        for(j in 1:series.count)
        {
            if(i >= j)
            {
                tmp.array[j,i] = params[iter]
                iter = iter + 1
            }
        }
    }
    buff.par[[1]] = tmp.array

    ## following loop initalizes the ARCH and GARCH parameter matrices respectively
    for(count in 1:(order[2] + order[1]))
    {
        buff.par[[count + 1]] = array(params[(count.triangular(series.count) + 1 + (count - 1) * series.count^2):(count.triangular(series.count) + 1 + series.count^2 + (count - 1) * series.count^2)], dim = c(series.count, series.count));
    }

    ## calculate the transposes of the parameter matrices
    buff.par.transposed = list()
    for(count in 1:length(buff.par))
    {
        buff.par.transposed[[count]] = t(buff.par[[count]])
    }

    ## compute the generalised kronecker product sums
    kronecker.sum = 0
    for(count in 2:total.par.matrices)
    {
        kronecker.sum = kronecker.sum + kronecker(buff.par[[count]], buff.par[[count]])
    }

    ## compute the eigenvalues
    tmp.svd = svd(kronecker.sum)
    eigenvalues = tmp.svd$d


    ## compute the unconditional covariance matrix
    numerat = t(buff.par[[1]]) %*% buff.par[[1]]
    dim(numerat) = c(series.count^2,1)
    denom = solve(diag(rep(1, series.count^2)) - kronecker.sum)
    sigma = denom %*% numerat
    dim(sigma) = c(series.count, series.count)

    ##  simulate a two-dimensional normal white noise process:
    T = T + 50    # the first 50 are later to be discarded
    nu = rnorm(series.count * T)
    nu = array(nu, dim = c(series.count, T))


    ## construct the simulated BEKK process
    HLAGS = list()
    for(count in 1:max(order))
    {
        HLAGS[[count]] = array(rep(0, series.count^2), dim = c(series.count,series.count))
        diag(HLAGS[[count]]) = 1
    }

    H = array(rep(0, series.count^2), dim = c(series.count, series.count))
    cor = numeric()
    eps.list = list()            # declare the estimated standard deviation series
    for(i in 1:series.count)
    {
        eps.list[[i]] = numeric()
    }

    sd = list()            # declare the estimated standard deviation series
    for(i in 1:series.count)
    {
        sd[[i]] = numeric()
    }


    ## initialize the first instances of the time series where
    ## HLAGS are not available
    for(count in 1:max(order))
    {
        for(i in 1:series.count)
        {
            eps.list[[i]][count] = 0
        }
    }

    eps = array(rep(0, series.count), dim = c(series.count,1))
    CTERM = buff.par.transposed[[1]] %*% buff.par[[1]] # the C'C term
    for(count in (max(order) + 1):T)
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
        H = CTERM
        ord1 = 1
        for(tmp.count in 1:(order[2] + order[1]))
        {
            if(tmp.count <= order[2])
            {
                ## ARCH EFFECT
                tmp.arr = numeric()
                for(scount in 1:series.count)
                {
                    tmp.arr[scount] = eps.list[[scount]][count - tmp.count]
                }
                eps = array(tmp.arr, dim = c(series.count,1))
                H = H + buff.par.transposed[[(tmp.count + 1)]] %*% eps %*% t(eps) %*% buff.par[[(tmp.count + 1)]]
            }
            else
            {
                ## GARCH EFFECT
                H = H + buff.par.transposed[[(tmp.count + 1)]] %*% HLAGS[[ord1]]  %*% buff.par[[(tmp.count + 1)]]
                ord1 = ord1 + 1
            }
        }

        svdH = svd(H)
        sqrtH = svdH$u %*% diag(sqrt(svdH$d)) %*% t(svdH$v)
        eps = sqrtH %*% nu[,count]

        ##cor[count] = H[1,2]/(sqrt(H[1,1] * H[2,2]))
        for(i in 1:series.count)
        {
            sd[[i]][count] = sqrt(H[i, i])
            eps.list[[i]][count] = eps[i,1]
        }
    }


    names(order) <- c("GARCH component", "ARCH component")
    names(buff.par) <- as.integer(seq(1, order[1] + order[2] + 1))

    nu = nu[51:T]

    for(i in 1:series.count)
    {
        sd[[i]] = sd[[i]][51:T]
        eps.list[[i]] = eps.list[[i]][51:T]
    }
    T = T - 50

    retval <- list(
        length = T,
        series.count = series.count,
        order = order,
        params = params,
        true.params = buff.par,
        eigenvalues = eigenvalues,
        uncond.cov.matrix = sigma,
        white.noise = nu,
        eps = eps.list,
        cor = cor,
        sd = sd
    )

    class(retval) <- "simulateBEKK"

    cat("Class attributes are accessible through following names:\n")
    cat(names(retval), "\n")

    return(retval)
}


##' Estimate MGARCH-BEKK processes
##'
##' Provides the MGARCH-BEKK estimation procedure.
##'
##' \code{BEKK} estimates a \code{BEKK(p,q)} model, where \code{p}
##' stands for the GARCH order, and \code{q} stands for the ARCH
##' order.
##'
##' @param eps Data frame holding time series.
##' @param order BEKK(p, q) order. An integer vector of length 2
##'     giving the orders of the model to be fitted. \code{order[2]}
##'     refers to the ARCH order and \code{order[1]} to the GARCH
##'     order.
##' @param params Initial parameters for the \code{optim} function.
##' @param fixed Vector of parameters to be fixed.
##' @param method The method that will be used by the \code{optim}
##'     function.
##' @param verbose Indicates if we need verbose output during the
##'     estimation.
##' @return Estimation results packaged as \code{BEKK} class
##'     instance. \describe{
##'         \item{eps}{a data frame contaning all time series}
##'         \item{length}{length of the series}
##'         \item{order}{order of the BEKK model fitted}
##'         \item{estimation.time}{time to complete the estimation process}
##'         \item{total.time}{time to complete the whole routine within the mvBEKK.est process}
##'         \item{estimation}{estimation object returned from the optimization process, using \code{optim}}
##'         \item{aic}{the AIC value of the fitted model}
##'         \item{est.params}{list of estimated parameter matrices}
##'         \item{asy.se.coef}{list of asymptotic theory estimates of standard errors of estimated parameters}
##'         \item{cor}{list of estimated conditional correlation series}
##'         \item{sd}{list of estimated conditional standard deviation series}
##'         \item{H.estimated}{list of estimated series of covariance matrices}
##'         \item{eigenvalues}{estimated eigenvalues for sum of Kronecker products}
##'         \item{uncond.cov.matrix}{estimated unconditional covariance matrix}
##'         \item{residuals}{list of estimated series of residuals}
##'     }
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
##' @examples
##' ## Simulate series:
##' simulated <- simulateBEKK(2, 1000, c(1,1))
##'
##' ## Prepare the matrix:
##' simulated <- do.call(cbind, simulated$eps)
##'
##' ## Estimate with default arguments:
##' estimated <- BEKK(simulated)
##'
##' \dontrun{
##' ## Show diagnostics:
##' diagnoseBEKK(estimated)
##' }
##'
##' @import mvtnorm
##' @import tseries
##' @import stats
##' @useDynLib mgarchBEKK
##' @export
BEKK <- function(eps, order  = c(1,1), params = NULL, fixed  = NULL, method = "BFGS", verbose = F) {
    ## TODO: Check the import statements in the preamble.
    ## TODO: Simplify count.triangular function.
    ## TODO: Use proper logger instead of out function.
    count.triangular <- function(dimension){
        if(dimension <= 0){
            0
        }
      else{
          dimension + count.triangular(dimension - 1)
      }
    }

    if(verbose == T){
        out <- function(...){
            cat(...)
        }
    }
    else{
        out <- function(...) { }
    }

    ## get the length and the number of the series
    series.length = length(eps[,1])
    series.count  = length(eps[1,])

    ## check the given order
    ## orders should be integers
    if(order[1] != as.integer(order[1]) || order[2] != as.integer(order[2]))
    {
        stop("order property should contain integer values")
    }

    ## GARCH effect could be set to 0, but, ARCH should be greater than 0
    if(order[1] < 0 || order[2] < 1)
    {
        stop("BEKK(",order[1],",",order[2],") is not implemented.")
    }

    ## construct the paramters list.
    ## first get the length of the parameter list
    params.length = count.triangular(series.count) + (order[2] * series.count^2) + (order[1] * series.count^2)

    if(is.null(params))
    {
        ## TODO
        ## these are meaningless parameters.
        ## set some useful initial parameters.
        params = c(1,0,1,0,0,1)
        params = c(params, rep(0.1, params.length - 6))
        out("\nWarning: initial values for the parameters are set to:\n\t", params,"\n")
    }
    else if(length(params) != params.length)
    {
        stop("Length of the initial parameter list doesn't conform required length. There should be ", params, " parameters in total")
    }

    ## check the given fixed parameters
    if(!is.null(fixed))
    {
        ## check the format of the fixed parameters
        if(
        (!is.array(fixed)) ||
        (dim(fixed)[1] != 2) ||
        (length(fixed[1,]) != length(fixed[2,])))
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
        for(i in 1:length(fixed[1,]))
        {
            fake.params[fixed[1,][i]] = NA
        }
        fake.params = na.omit(fake.params)
    }

    ## parameters seem appropriate
    ## define the loglikelihood function
    loglikelihood.C <- function(params)
    {
        loglikelihood.C <- .C("loglikelihood",
                              as.vector(params, mode = "double"),
                              as.vector(fixed[1,], mode = "integer"),
                              as.vector(fixed[2,], mode = "double"),
                              as.integer(length(fixed[1,])),
                              as.vector(t(eps)), # funny: transpose the time series
                              as.integer(series.count),
                              as.integer(series.length),
                              as.vector(order, mode = "integer"),
                              retval = 0.0,
                              PACKAGE = "mgarchBEKK"
                              )

        if(is.nan(loglikelihood.C$retval) == T)
        {
            nonusedret = 1e+100

        }
        else
        {
            nonusedret = loglikelihood.C$retval
        }
        nonusedret
    }

    ## begin estimation process

    ## first log the start time
    start = Sys.time()
    out("* Starting estimation process.\n")
    out("* Optimization Method: '", method, "'\n")

    ## call the optim function
    estimation = optim(fake.params, loglikelihood.C, method = method, hessian = T)

    ## estimation completed
    out("* Estimation process completed.\n")

    ## log estimation time
    est.time = difftime(Sys.time(), start)

    ## calculate the AIC
    ## it is estimation value + number of estimated parameters (punishment :))
    aic = estimation$value + (length(params) - length(fixed[1,]))

    ## following script will prepare an object that holds the estimated
    ## parameters and some useful diagnostics data like estimated correlation,
    ## standard deviation, eigenvalues and so on.

    ## TODO
    ## estimation$hessian is non-existing if fixed parameter list contains all the
    ## paramters to be estimated. That is that the estimation procedure gets no parameters,
    ## thus, there is no errors... Fix it... How?
    ## Whether encapsulate with an "if" statement, probably not efficient,
    ## or give a fake hessian
    ##
    ## give a fake hessian
    if(length(fake.params) == 0)
    {
        estimation$hessian = matrix(rep(0.1, series.count^2), nrow = series.count, ncol = series.count)
    }

    ## get the hessian matrix and grap the diagonal
    inv.hessian.mat = solve(estimation$hessian)

    diag.inv.hessian = sqrt(abs(diag(inv.hessian.mat)))
    if(length(which(diag(inv.hessian.mat) < 0)) == 0)
    {
        warning("negative inverted hessian matrix element")
    }

    ## fix the asymptotic-theory standard errors of the
    ## coefficient estimates with fixed parameters
    if(!is.null(fixed))
    {
        temp.diag.inv.hessian = numeric()
        shifted = 0
        for(count in 1:params.length)
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
    parnum = 1 + order[1] + order[2]    # calculate number of paramater matrices
    asy.se.coef = list()            # declare the asy.se.coef matrices list

    ############################################################################
    ## CRITICAL!!!
    ## since we are not bidimensional anymore, be careful!!!
    ############################################################################
    ## first initialize the first asy.se.coef matrix, corresponding to the C matrix
    tmp.array = array(rep(0, series.count^2), dim = c(series.count, series.count))
    tmp.array[!lower.tri(tmp.array)] = diag.inv.hessian[1:length(which(!lower.tri(tmp.array) == T))]
    asy.se.coef[[1]] = tmp.array

    ## following loop initalizes the ARCH and GARCH parameter matrices respectively
    for(count in 1:(parnum - 1))
    {
        ## !! a bit hard to follow
        asy.se.coef[[count + 1]] = array(diag.inv.hessian[(count.triangular(series.count) + 1 + (count - 1) * series.count^2):(count.triangular(series.count) + 1 + series.count^2 + (count - 1) * series.count^2)], dim = c(series.count, series.count));
    }
    #
    buff.par = list()        # declare the parameter list

    ## shift the fixed parameters inside the estimated parameters
    if(!is.null(fixed))
    {
        estim.params = numeric()
        shifted = 0
        for(count in 1:params.length)
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
    tmp.array = array(rep(0, series.count^2), dim = c(series.count, series.count))
    tmp.array[!lower.tri(tmp.array)] = estim.params[1:length(which(!lower.tri(tmp.array) == T))]

    buff.par[[1]] = tmp.array

    ## following loop initalizes the ARCH and GARCH parameter matrices respectively
    for(count in 1:(parnum - 1))
    {
        ## !! a bit hard to follow
        buff.par[[count + 1]] = array(estim.params[(count.triangular(series.count) + 1 + (count - 1) * series.count^2):(count.triangular(series.count) + 1 + series.count^2 + (count - 1) * series.count^2)], dim = c(series.count, series.count));
    }

    ## calculate the transposes of the parameter matrices
    buff.par.transposed = lapply(buff.par, t)

    ## start diagnostics
    out("* Starting diagnostics...\n")
    out("* Calculating estimated:\n")
    out("*\t1. residuals,\n")
    out("*\t2. correlations,\n")
    out("*\t3. standard deviations,\n")
    out("*\t4. eigenvalues.\n")

    HLAGS = list()        # list of H lags that will be used later in the MGARCH implementation
    for(count in 1:order[1])
    {
        ## TODO:check intial values (currently 1's on the diagonal)
        HLAGS[[count]] = array(rep(0, series.count^2), dim = c(series.count,series.count))
        diag(HLAGS[[count]]) = 1
    }

    residuals = list()
    for(i in 1:series.count)
    {
        residuals[[i]] = numeric()
    }

    ## initialize the first residuals we are not able to calculate
    for(count in 1:max(order))
    {
        for(i in 1:series.count)
        {
            residuals[[i]][count] = 0
        }
    }

    resid = array(rep(0,series.count), dim = c(series.count,1)) # declare a temporary residuals buffer

    ## calculate eigenvalues
    ## TODO:
    ## Angi says that following is not true according to Bauwens, Laurent, Rombouts Paper.
    temp = 0
    for(count in 2:parnum)
    {
        temp = temp + kronecker(buff.par[[count]], buff.par[[count]])
    }
    eigenvalues = svd(temp)$d

    ##################################################################
    ### TODO:
    ### FROM NOW ON, HELP NEEDED
    ### ASK HARALD HOCA!
    ##################################################################

    ## compute the unconditional covariance matrix
    numerat = t(buff.par[[1]]) %*% buff.par[[1]]
    dim(numerat) = c(series.count^2,1)
    denom = solve(diag(rep(1, series.count^2)) - temp)
    sigma = denom %*% numerat
    dim(sigma) = c(series.count, series.count)

    H = cov(eps)    # to initialize, use the covariance matrix of the series
    H.estimated = lapply(1:series.length, function(x){H})

    cor = list()        # declare the estimated correlation series
    for(i in 1:series.count)
    {
        cor[[i]] = list()
        for(j in 1:series.count)
        {
            cor[[i]][[j]] = numeric()
        }
    }
    sd = list()        # declare the estimated standard deviation series
    for(i in 1:series.count)
    {
        sd[[i]] = numeric()
    }

    eps.est = array(rep(0,series.count), dim = c(series.count,1))    # declare a temporary eps buffer

    CTERM = buff.par.transposed[[1]] %*% buff.par[[1]] # calculate the C'C term

    out("* Entering Loop...");
    for(count in (max(order) + 1):series.length) # cruical loop! initializing diagnostics data
    {
        ## do the swap calculation for H terms
        if(order[1] >= 2)
        {
            for(tmp.count in order[1]:2)
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
                H = H + buff.par.transposed[[tmp.count + 1]] %*% as.matrix(eps[count - tmp.count,]) %*% as.matrix(t(eps[count - tmp.count,])) %*% buff.par[[tmp.count + 1]]
            }
            else
            {
                ## GARCH EFFECT (G TERM)
                H = H + buff.par.transposed[[tmp.count + 1]] %*% HLAGS[[ord1]]  %*% buff.par[[tmp.count + 1]]
                ord1 = ord1 + 1
            }
        }

        ## TODO add appropriate comments for following assignments and calculations
        H.estimated[[count]] = H
        svdH = svd(H)
        sqrtH = svdH$u %*% diag(sqrt(svdH$d)) %*% t(svdH$v)

        invsqrtH = solve(sqrtH)
        resid = invsqrtH %*% as.matrix(eps[count,])
        for(i in 1:series.count)
        {
            residuals[[i]][count] = resid[i,1]
        }

        ## TODO: check
        for(i in 1:series.count)
        {
            for(j in 1:series.count)
            {
                cor[[i]][[j]][count] = H[i,j] / sqrt(H[i,i] * H[j,j])
            }
        }
        for(i in 1:series.count)
        {
            sd[[i]][count] = sqrt(H[i,i])
        }
    }

    ## diagnostics ready
    out("Diagnostics ended...\n")

    names(order) <- c("GARCH component", "ARCH component")
    names(buff.par) <- as.integer(seq(1, parnum))

    retval <- list(
        eps = eps,
        series.length = series.length,
        estimation.time = est.time,
        total.time = difftime(Sys.time(), start),
        order = order,
        estimation = estimation,
        aic = aic,
        asy.se.coef = asy.se.coef,
        est.params = buff.par,
        cor = cor,
        sd = sd,
        H.estimated = H.estimated,
        eigenvalues = eigenvalues,
        uncond.cov.matrix = sigma,
        residuals = residuals
    )

    class(retval) = "BEKK"

    out("Class attributes are accessible through following names:\n")
    out(names(retval), "\n")

    return(retval)
}


##' Diagnose BEKK process estimation
##'
##' Provides diagnostics for a BEKK process estimation.
##'
##' This procedure provides console output and browsable plots for a
##' given BEKK process estimation. Therefore, it is meant to be
##' interactive as the user needs to proceed by pressing \code{c} on
##' the keyboard to see each plot one-by-one.
##'
##' @param estimation The return value of the \code{mvBEKK.est} function
##' @return Nothing special
##'
##' @examples
##' ## Simulate series:
##' simulated = simulateBEKK(2, 1000, c(1,1))
##'
##' ## Prepare the matrix:
##' simulated = do.call(cbind, simulated$eps)
##'
##' ## Estimate with default arguments:
##' estimated = BEKK(simulated)
##'
##' \dontrun{
##' ## Show diagnostics:
##' diagnoseBEKK(estimated)
##' }
##'
##' @import stats
##' @import graphics
##' @import grDevices
##' @export
diagnoseBEKK <- function(estimation)
{
    cat("\tNumber of estimated series : ", length(estimation$eps),     "\n")
    cat("\tLength of estimated series : ", estimation$series.length,   "\n")
    cat("\tEstimation Time            : ", estimation$estimation.time, "\n")
    cat("\tTotal Time                 : ", estimation$total.time,      "\n")
    cat("\tBEKK order                 : ", estimation$order,           "\n")
    cat("\tEigenvalues                : ", estimation$eigenvalues,     "\n")
    cat("\taic                        : ", estimation$aic,             "\n")
    cat("\tunconditional cov. matrix  : ", estimation$uncond.cov.mat,  "\n")

    for(i in 1:length(estimation$eps[1,]))
    {
        cat("\tvar(resid", i, ")                : ", var(estimation$residuals[[i]]),      "\n")
        cat("\tmean(resid", i, ")               : ", mean(estimation$residuals[[i]]),     "\n")
    }

    cat("\tEstimated parameters       :\n\n")
    cat("\tC estimates:\n")
    print(estimation$est.params[[1]])

    if(estimation$order[2] > 0)
    {
        cat("\n\tARCH estimates:\n")
        for(count in 1:estimation$order[2])
        {
            print(estimation$est.params[[count + 1]])
        }
    }
    else
    {
        count = 0
    }

    if(estimation$order[1] > 0)
    {
        cat("\n\tGARCH estimates:\n")
        for(count2 in 1:estimation$order[1])
        {
            print(estimation$est.params[[(count + 1) + count2]])
        }
    }

    cat("\n\tasy.se.coef                : \n\n")
    cat("\tC estimates, standard errors:\n")
    print(estimation$asy.se.coef[[1]])

    if(estimation$order[2] > 0)
    {
        cat("\n\tARCH estimates, standard errors:\n")
        for(count in 1:estimation$order[2])
        {
            print(estimation$asy.se.coef[[count + 1]])
        }
    }
    else
    {
        count = 0
    }

    if(estimation$order[1] > 0)
    {
        cat("\n\tGARCH estimates, standard errors:\n")
        for(count2 in 1:estimation$order[1])
        {
            print(estimation$asy.se.coef[[(count + 1) + count2]])
        }
    }

    ##    plot(
    ##            min(min(estimation$resid1),min(estimation$resid2)):max(max(estimation$resid1),max(estimation$resid2)),
    ##            min(min(estimation$resid1),min(estimation$resid2)):max(max(estimation$resid1),max(estimation$resid2)),
    ##            type = "n",
    ##            xlab = "resid1",
    ##            ylab = "resid2"
    ##        )
    ##
    ##    points(estimation$resid1, estimation$resid2, pch = 21)

    for(i in 1:length(estimation$eps[1,]))
    {
        plot(estimation$residuals[[i]])
        browser()
        dev.off()
    }
}
