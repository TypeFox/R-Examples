#' Testing the equality of critical points
#'@description This function can be used to test the equality of the
#' \eqn{M} critical points estimated from the respective level-specific curves.
#' @param formula An object of class \code{formula}: a sympbolic 
#' description of the model to be fitted. The details of model 
#' specification are given under 'Details'.
#'@param data A data frame or matrix containing the model response variable
#' and covariates required by the \code{formula}.
#' @param der Number which determines any inference process. 
#' By default \code{der} is \code{NULL}. If this term is \code{0}, 
#' the testing procedures is applied for the estimate. If it is \code{1} or
#' \code{2}, it is designed for the first or second derivative, respectively.
#' @param smooth Type smoother used: \code{smooth = "kernel"} for local polynomial
#' kernel smoothers and \code{smooth = "splines"} for splines using the 
#' \code{mgcv} package.
#' @param weights Prior weights on the data.
#' @param nboot Number of bootstrap repeats.
#' @param h0 The kernel bandwidth smoothing parameter for the global effect (see
#' references for more details at the estimation). Large values of the bandwidth lead
#' to smoothed estimates; smaller values of the bandwidth lead lo undersmoothed estimates. 
#' By default, cross validation is used to obtain the bandwidth.
#' @param h The kernel bandwidth smoothing parameter for the partial effects.
#' @param nh Integer number of equally-spaced bandwidth on which the
#' \code{h} is discretised, to speed up computation.
#' @param kernel A character string specifying the desired kernel. 
#' Defaults to \code{kernel = "epanech"}, where the Epanechnikov
#' density function kernel will be used. Also, several types of kernel funcitons 
#' can be used:  triangular and Gaussian density function, 
#' with \code{"triang"} and \code{"gaussian"} term, respectively.
#' @param p Degree of polynomial to be used. Its value must be the value of
#' derivative + 1. The default value is 3 due to the function
#' returns the estimation, first and second derivative.
#' @param kbin Number of binning nodes over which the function 
#' is to be estimated.
#' @param rankl Number or vector specifying the minimum value for the
#' interval at which to search the \code{x} value which maximizes the
#' estimate, first or second derivative  (for each level). The default
#' is the minimum data value.
#' @param ranku Number or vector specifying the maximum value for the
#' interval at which to search the \code{x} value which maximizes the
#' estimate, first or second derivative  (for each level). The default
#' is the maximum data value.
#' @param seed Seed to be used in the bootstrap procedure.
#' @param cluster A logical value. If  \code{TRUE} (default), the
#'  bootstrap procedure is  parallelized (only for \code{smooth = "splines"}.
#'   Note that there are cases 
#'  (e.g., a low number of bootstrap repetitions) that R will gain in
#'  performance through serial computation. R takes time to distribute tasks
#'  across the processors also it will need time for binding them all together
#'  later on. Therefore, if the time for distributing and gathering pieces
#'  together is greater than the time need for single-thread computing, it does
#'  not worth parallelize.
#'@param ncores An integer value specifying the number of cores to be used
#' in the parallelized procedure. If \code{NULL} (default), the number of cores 
#' to be used is equal to the number of cores of the machine - 1.
#' @param \ldots Other options.
#' 
#' 
#' @details \code{localtest} can be used to test the equality of the 
#' \eqn{M} critical points estimated from the respective level-specific curves. 
#' Note that, even if the curves and/or their derivatives are different, it is 
#' possible for these points to be equal. 
#' 
#' For instance, taking the maxima of the first derivatives into account, 
#' interest lies in testing the following null hypothesis
#' 
#' \deqn{H_0: x_{01} = \ldots = x_{0M}}
#' 
#' versus the general  alternative  
#' 
#' \deqn{H_1: x_{0i} \ne x_{0j}  \quad {\rm{for}} \quad {\rm{some}} \quad 
#' \emph{i}, \emph{j} \in \{ 1, \ldots, M\}.}
#' 
#' The above hypothesis is true if \eqn{d=x_{0j}-x_{0k}=0} where 
#' \deqn{ (j,k)= argmax \quad (l,m) \quad \{1 \leq l<m \leq M\} \quad |x_{0l}-x_{0m}|, }
#' 
#' otherwise  \eqn{H_0} is false. It is important to highlight that, in practice,
#' the true \eqn{x_{0j}} are not known, and consequently neither is \eqn{d}, 
#' so an estimate \eqn{\hat d = \hat x_{0j}-\hat x_{0k}} is used, where, 
#' in general, \eqn{\hat x_{0l}} are the estimates of \eqn{x_{0l}} based on the 
#' estimated curves \eqn{\hat m_l} with \eqn{l = 1, \ldots , M}. 
#' 
#' Needless to say, 
#' since \eqn{\hat d} is only an estimate of the true \eqn{d}, the sampling 
#' uncertainty of these estimates needs to be acknowledged. Hence, a confidence 
#' interval \eqn{(a,b)} is created for \eqn{d} for a specific level of 
#' confidence (95\%).  Based on this, the null hypothesis is rejected if  
#' zero is not contained in the interval.
#' 
#' Note that if this hypothesis is rejected (and the factor has more than 
#' two levels), one option could be to use the \code{maxp.diff} function in 
#' order to obtain the differences between each pair of factor's levels.
#' 
#' Note that the models fitted by \code{localtest} function are specified 
#' in a compact symbolic form. The \~ operator is basic in the formation 
#' of such models. An expression of the form \code{y ~ model}  is interpreted as 
#' a specification that the response \code{y} is modelled by a predictor 
#' specified symbolically by \code{model}. The possible terms consist of a 
#' variable name or a variable name and a factor name separated by : operator. 
#' Such a term is interpreted as the interaction of the continuous variable and 
#' the factor. However, if \code{smooth = "splines"}, the formula is based on the function
#' formula.gam of the mgcv package.

#'@return The estimate of \eqn{d} value is returned and its confidence interval 
#'for a specific-level of confidence, i.e. 95\%. Additionally, it is shown 
#'the decision, accepted or rejected,  of the local test. Based on the null 
#'hypothesis is rejected if a zero value is not within the interval. 
#'
#'@author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#'
#' @references 
#' Sestelo, M. (2013). Development and computational implementation of 
#' estimation and inference methods in flexible regression models. 
#' Applications in Biology, Engineering and Environment. PhD Thesis, Department
#' of Statistics and O.R. University of Vigo.
#' 
#'@examples
#' library(npregfast)
#' data(barnacle)
#' localtest(DW ~ RC : F, data = barnacle, der = 1, seed = 130853, nboot = 100)
#' 
#' @useDynLib npregfast localtest_
#' @importFrom stats na.omit runif
#' @importFrom mgcv interpret.gam gam predict.gam
#' @importFrom sfsmisc D1D2
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores
#' @importFrom foreach foreach %dopar%
#' @export




localtest <- function(formula, data = data, der, smooth = "kernel", weights = NULL, 
                      nboot = 500, h0 = -1.0, h = -1.0, nh = 30, kernel = "epanech", 
                      p = 3, kbin = 100, rankl = NULL, ranku = NULL, seed = NULL,
                      cluster = TRUE, ncores = NULL, ...) {
  
  if(kernel == "gaussian")  kernel <- 3
  if(kernel == "epanech")   kernel <- 1
  if(kernel == "triang")    kernel <- 2
  
  if (missing(der)) {
    stop("Argument \"der\" is missing, with no default")
  }
  
  if (missing(formula)) {
    stop("Argument \"formula\" is missing, with no default")
  }
  if (missing(data)) {
    stop("Argument \"data\" is missing, with no default")
  }
  
  if(!isTRUE(der %in% c(0, 1, 2))) {
    stop("",paste(der)," is not a r-th derivative implemented, only 
         permitted 0, 1 or 2.")
  }
  
  if (!(kernel %in% 1:3)) {
    stop("Kernel not suported")
  }
  
  if (!(smooth %in% c("kernel", "splines"))) {
    stop("Smoother not suported")
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (isTRUE(cluster) & smooth == "splines") {
    if (is.null(ncores)) {
      num_cores <- detectCores() - 1
    }else{
      num_cores <- ncores
    }
    registerDoParallel(cores = num_cores)
  }
  
  ncmax <- 5
  c2 <- NULL
  # if(is.null(seed)) seed <- -1
  
  
  
  if (smooth != "splines") {
    
    ffr <- interpret.frfastformula(formula, method = "frfast")
    varnames <- ffr$II[2, ]
    aux <- unlist(strsplit(varnames,split = ":"))
    varnames <- aux[1]
    if (unlist(strsplit(varnames,split = ""))[1] == "s") {
      stop("Argument \"formula\" is wrong specified, see details of
           model specification in 'Details' of the frfast help." )
    }
    namef <- aux[2]
    if (length(aux) == 1) {f <- NULL}else{f <- data[ ,namef]}
    newdata <- data
    data <- na.omit(data[ ,c(ffr$response, varnames)])
    #newdata <- na.omit(newdata[ ,varnames])
    n <- nrow(data)
    
  }else{
    ffr <- interpret.gam(formula)
    varnames <- ffr$pred.names[1]
    if (":" %in% unlist(strsplit(ffr$fake.names,split = ""))) {
      stop("Argument \"formula\" is wrong specified, see details of
             model specification in 'Details' of the frfast help." )
    }
    
    namef <- ffr$pred.names[2]
    if (length(ffr$pred.names) == 1) {f <- NULL}else{f <- data[ ,namef]}
    newdata <- data
    if (length(ffr$pred.names) == 1) {
      data <- na.omit(data[ ,c(ffr$response, varnames)])
    }else{
      data <- na.omit(data[ ,c(ffr$response, varnames, namef)])
    }
    #newdata <- na.omit(newdata[ ,varnames])
    n <- nrow(data)
  }
  
  
  if (is.null(f)) f <- rep(1, n)
  etiquetas <- unique(f)
  nf <- length(etiquetas)
  
  
  if(nf == 1) {
    stop("Function not supported.
         There is not factor in the model.")
  }
  
  if(is.null(h0)){
    h0 <- -1.0
  }
  if(is.null(h)){
    h <- rep(-1.0, nf)
  }else{
    if(length(h) == 1) h <- rep(h, nf)
  }
  
  # Interesaria meter para las derivadas?
  
  if (is.null(weights)) {
    weights <- rep(1, n)
  } else {
    if (sum(weights) <= 0 || any(weights) < 0 || length(weights) != n) 
      stop("The specified weights are not correct")
  }
  
  
  if(is.null(c2)) c2 <- matrix(as.double(-1.0), ncmax, nf) 
  if(is.null(rankl)){
    rankl <- na.omit(as.vector(tapply(data[ ,varnames], f, min)))
  }else{
    if(length(rankl) == 1) rankl <- rep(rankl, nf)
  }
  if(is.null(ranku)){
    ranku <- na.omit(as.vector(tapply(data[ ,varnames], f, max)))
  }else{
    if(length(ranku) == 1) ranku <- rep(ranku, nf)
  } 
  
  
  
  if (smooth != "splines") {
    
    
    umatrix <- matrix(runif(n*nboot), ncol = nboot, nrow = n)
    
    
    localtest  <-.Fortran("localtest_",
                          f = as.integer(f),
                          x = as.double(data[,varnames]),
                          y = as.double(data[,ffr$response]),
                          w = as.double(weights),
                          n = as.integer(n),
                          h0 = as.double(h0),
                          h = as.double(h),
                          nh = as.integer(nh),
                          p = as.integer(p),
                          kbin = as.integer(kbin),
                          #fact = as.integer(c(1:nf)),
                          fact = unique(as.integer(f)),
                          #fact   =as.integer(c(1:nf))
                          nf = as.integer(nf),
                          kernel = as.integer(kernel),
                          nboot = as.integer(nboot),
                          pcmax = as.double(ranku), # rango de busqueda minimo
                          pcmin = as.double(rankl), # rango de busqueda maximo
                          r = as.integer(der),
                          D = as.double(rep(-1.0,1)),
                          Ci = as.double(rep(-1.0,1)),
                          Cs = as.double(rep(-1.0,1)),
                          seed = as.integer(seed),
                          umatrix = as.double(umatrix)
    )
    
    
    
    if (localtest$Ci <= 0 & 0 <= localtest$Cs) {
      decision <- "Acepted"
    } else {
      decision <- "Rejected"
    }
    res <- cbind(d = round(localtest$D, digits = 4), Lwr = round(localtest$Ci, digits = 4), 
                 Upr = round(localtest$Cs, digits = 4), Decision = decision)
    # class(res) <- 'localtest'
    
  }else{
    
    mainfun_localtest <- function(formula, data, weights){
      
      # grid
      xgrid <- seq(min(data[ ,varnames]), max(data[ ,varnames]), length.out = kbin)
      newd <- expand.grid(xgrid, unique(f))
      names(newd) <- ffr$pred.names
      
      # estimations
      p <- array(NA, dim = c(kbin, 3, nf))
      m <- gam(formula, weights = weights, data = data.frame(data, weights), ...)
      muhat <- as.vector(predict(m, newdata = newd, type = "response"))
      p[, 1, 1:nf] <- muhat
      #d1 <- apply(p, 3, function(z){D1ss(x = xgrid, y = z[, 1])})
      d1 <- apply(p, 3, function(z){D1D2(x = xgrid, y = z[, 1], deriv = 1)$D1})
      p[, 2, 1:nf] <- as.vector(d1)
      #d2 <- apply(p, 3, function(z){D2ss(x = xgrid, y = z[, 1])$y})
      d2 <- apply(p, 3, function(z){D1D2(x = xgrid, y = z[, 1], deriv = 2)$D2})
      p[, 3, 1:nf] <- as.vector(d2)
      
      
      # gridfino
      kfino <- 100
      xgridfino <- vapply(c(1:nf),
                          FUN = function(x){seq(rankl[x], ranku[x], length.out = kfino)}, 
                          FUN.VALUE = numeric(kfino))
      
      newdfino <- data.frame(as.vector(xgridfino), rep(unique(f), each = kfino))
      names(newdfino) <- ffr$pred.names
      
      # max
      muhatfino <- as.vector(predict(m, newdata = newdfino, type = "response"))
      
      pfino <- array(NA, dim = c(kfino, 3, nf))
      pfino[, 1, 1:nf] <- muhatfino
      
      aux <- data.frame(muhatfino, newdfino)
      # d1 <- by(aux, aux[, 3], function(z){D1ss(x = z[, 2], y = z[, 1])})
      d1 <- by(aux, aux[, 3], function(z){D1D2(x = z[, 2], y = z[, 1], deriv = 1)$D1})
      pfino[, 2, 1:nf] <- unlist(d1)
      #d2 <- by(aux, aux[, 3], function(z){D2ss(x = z[, 2], y = z[, 1])$y})
      d2 <- by(aux, aux[, 3], function(z){D1D2(x = z[, 2], y = z[, 1], deriv = 2)$D2})
      pfino[, 3, 1:nf] <- unlist(d2)
      
      
      iimax <- apply(pfino, 3:2, which.max)
      iicero <- apply(abs(pfino), 3:2, which.min)
      max <- matrix(NA, ncol = nf, nrow = 3)
      for (j in 1:nf) {
        max[1:2 , j] <- xgridfino[ t(iimax)[1:2, j], j]
        max[3 , j] <- xgridfino[ t(iicero)[3, j], j]
      }
      
      xmin <- min(max[der + 1, ])
      posmin <- which.min(max[der + 1, ])
      
      xmax <- max(max[der + 1, ])
      posmax <- which.max(max[der + 1, ])
      
      if (posmin < posmax){
        d <-  xmin - xmax
      }else{
        d <-  xmax - xmin
      }
      
      return(d)
    }
    
    
    
    d <- mainfun_localtest(formula, data = data, weights = weights)
    
    # bootstrap
    m <- gam(formula, weights = weights, data = data.frame(data, weights), ...)
    muhat <- as.vector(predict(m, type = "response"))
    err <- data[, ffr$response] - muhat
    err <- err - mean(err)
    yboot <- replicate(nboot, muhat + err *
                         sample(c(-sqrt(5) + 1, sqrt(5) + 1)/2, size = n,
                                replace = TRUE,
                                prob = c(sqrt(5) + 1, sqrt(5) - 1)/(2 * sqrt(5))))
    
    i <- NULL
    d_allboot <- foreach(i = 1:nboot) %dopar% {
      datab <- data
      datab[, ffr$response] <- yboot[, i]
      aux <- mainfun_localtest(formula, data = data.frame(datab, weights), 
                               weights = weights, ...)
      return(aux)
    }
    
    
    ci <- quantile(unlist(d_allboot), probs = c(0.025, 0.975), na.rm = TRUE)
    
    if (ci[1] <= 0 & 0 <= ci[2]) {
      decision <- "Acepted"
    } else {
      decision <- "Rejected"
    }
    res <- cbind(d = round(d, digits = 4), Lwr = round(ci[1], digits = 4), 
                 Upr = round(ci[2], digits = 4), Decision = decision)
    
    rownames(res) <- NULL
    
  }
  
  return(as.data.frame(res))
  
} 
