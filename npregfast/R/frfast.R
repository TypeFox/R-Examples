#' Fitting nonparametric models
#' 
#' @description This function is used to fit nonparametric models by
#' using local polynomial kernel smoothers or splines. These models can 
#' include or not  
#' factor-by-curve interactions. Additionally, a parametric 
#' model (allometric model) can be estimated (or not).
#' @param formula An object of class \code{formula}: a sympbolic 
#' description of the model to be fitted. The details of model 
#' specification are given under 'Details'.
#' @param data A data frame or matrix containing the model response
#' variable and covariates required by the \code{formula}.
#' @param model Type model used: \code{model = "np"}  nonparametric
#' regression model, 
#' \code{model = "allo"} the  allometric model.
#' @param smooth Type smoother used: \code{smooth = "kernel"} for local polynomial
#' kernel smoothers and \code{smooth = "splines"} for splines using the 
#' \code{mgcv} package.
#' @param h0 The kernel bandwidth smoothing parameter for the global effect (see
#' references for more details at the estimation). Large values of the bandwidth lead
#' to smoothed estimates; smaller values of the bandwidth lead lo undersmoothed estimates. 
#' By default, cross validation is used to obtain the bandwidth.
#' @param h The kernel bandwidth smoothing parameter for the partial effects. 
#' @param nh Integer number of equally-spaced bandwidth in which the
#' \code{h} is discretised, to speed up computation in the kernel-based regression.
#' @param weights Prior weights on the data.
#' @param kernel A character string specifying the desired kernel. 
#' Defaults to \code{kernel = "epanech"}, where the Epanechnikov
#' density function kernel will be used. Also, several types of kernel functons 
#' can be used:  triangular and Gaussian density function, 
#' with \code{"triang"} and \code{"gaussian"} term, respectively.
#' @param p Polynomial degree to be used in the kernel-based regression. Its 
#' value must be the value of
#' derivative + 1. The default value is 3, returning 
#' the estimation, first and second derivative.
#' @param kbin Number of binning nodes over which the function 
#' is to be estimated.
#' @param nboot Number of bootstrap repeats. Defaults to 500 bootstrap repeats. 
#' The wild bootstrap is used when \code{model = "np"} and the simple bootstrap 
#' when \code{model = "allo"}.
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
#'  bootstrap procedure is  parallelized (only for \code{smooth = "splines"}).
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
#'  
#' @details The models fitted by \code{frfast} function are specified 
#' in a compact symbolic form. The \~ operator is basic in the formation 
#' of such models. An expression of the form \code{y ~ model}  is interpreted as 
#' a specification that the response \code{y} is modelled by a predictor 
#' specified symbolically by \code{model}. The possible terms consist of a 
#' variable name or a variable name and a factor name separated by : operator. 
#' Such a term is interpreted as the interaction of the continuous variable and 
#' the factor. However, if \code{smooth = "splines"}, the formula is based on the function
#' formula.gam of the mgcv package.
#' @return An object is returned with the following elements:
#' \item{x}{Vector of values of the grid points at which model is to 
#' be estimate.}
#' \item{p}{Matrix of values of the grid points at which to compute the 
#' estimate, their first and second derivative.}
#' \item{pl}{Lower values of  95\% confidence interval for the estimate, 
#' their first and second derivative.}
#' \item{pu}{Upper values of  95\% confidence interval for the estimate, 
#' their first and second derivative.}
#' \item{diff}{Differences between the estimation values of a couple of 
#' levels (i. e. level 2 - level 1). The same procedure for their first
#' and second derivative.}
#' \item{diffl}{Lower values of 95\% confidence interval for the differences 
#' between the estimation values of a couple of levels. It is performed 
#' for their first and second derivative.}
#' \item{diffu}{Upper values of 95\% confidence interval for the differences 
#' between the estimation values of a couple of levels. It is performed for 
#' their first and second derivative.}
#' \item{nboot}{Number of bootstrap repeats.}
#' \item{n}{Sample size.}
#' \item{dp}{Degree of polynomial to be used.}
#' \item{h0}{The kernel bandwidth smoothing parameter for the global effect.}
#' \item{h}{The kernel bandwidth smoothing parameter for the partial effects.}
#' \item{fmod}{Factor's level for each data.}
#' \item{xdata}{Original x values.}
#' \item{ydata}{Original y values.}
#' \item{w}{Weights on the data.}
#' \item{kbin}{Number of binning nodes over which the function is to 
#' be estimated.}
#' \item{nf}{Number of levels.}
#' \item{max}{Value of covariate \code{x} which maximizes the  estimate, 
#' first or second derivative.}
#' \item{maxu}{Upper value of 95\% confidence interval for the 
#' value \code{max}.}
#' \item{maxl}{Lower value of 95\% confidence interval for the 
#' value \code{max}.}
#' \item{diffmax}{Differences between the estimation of \code{max} for a 
#' couple of levels (i. e. level 2 - level 1). The same procedure for their 
#' first and second derivative.}
#' \item{diffmaxu}{Upper value of 95\% confidence interval for the value 
#' \code{diffmax}.}
#' \item{diffmaxl}{Lower value of 95\% confidence interval for the value 
#' \code{diffmax}.}
#' \item{repboot}{Matrix of values of the grid points at which to compute 
#' the estimate, their first and second derivative for each bootstrap repeat.}
#' \item{rankl}{Maximum value for the interval at which to search the 
#' \code{x} value which maximizes the estimate, first or second derivative  
#' (for each level). The default is the maximum data value.}
#' \item{ranku}{Minimum value for the interval at which to search the 
#' \code{x} value which maximizes the estimate, first or second derivative  
#' (for each level). The default is the minimum data value.}
#' \item{nmodel}{Type model used: \code{nmodel = 1} the nonparametric model, 
#' \code{nmodel = 2} the allometric model.}
#' \item{label}{Labels of the variables in the model.}
#' \item{numlabel}{Number of labels.}
#' \item{kernel}{A character specifying the derised kernel.}
#' \item{a}{Estimated coefficient in the case of fitting an allometric model.}
#' \item{al}{Lower value of 95\% confidence interval for the value of \code{a}.}
#' \item{au}{Upper value of 95\% confidence interval for the value of \code{a}.}
#' \item{b}{Estimated coefficient in the case of fitting an allometric model.}
#' \item{bl}{Lower value of 95\% confidence interval for the value of \code{b}.}
#' \item{bu}{Upper value of 95\% confidence interval for the value of \code{b}.}
#' \item{name}{Name of the variables in the model.}
#' \item{formula}{A sympbolic description of the model to be fitted.}
#' \item{nh}{Integer number of equally-spaced bandwidth on which the
#' \code{h} is discretised.}
#' \item{r2}{Coefficient of determination (in the case of the allometric model).}
#' \item{smooth}{Type smoother used.}
#' \item{cluster}{Is the procedure parallelized? (for splines smoothers).}
#' \item{ncores}{Number of cores used in the parallelized procedure? (for splines smoothers).}
#' 
#' 
#' @author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#' 
#' @references 
#' Sestelo, M. (2013). Development and computational implementation of 
#' estimation and inference methods in flexible regression models. 
#' Applications in Biology, Engineering and Environment. PhD Thesis, Department
#' of Statistics and O.R. University of Vigo.
#' 
#' @examples
#' library(npregfast)
#' data(barnacle)
#' 
#' # Nonparametric regression without interactions
#' fit <- frfast(DW ~ RC, data = barnacle, nboot = 100, smooth = "kernel") 
#' fit
#' summary(fit)
#' 
#' # using  splines
#' #fit <- frfast(DW ~ s(RC), data = barnacle, nboot = 100, 
#' #smooth = "splines", cluster = TRUE, ncores = 2) 
#' #fit
#' #summary(fit)
#' 
#' 
#' # Change the number of binning nodes and bootstrap replicates
#' fit <- frfast(DW ~ RC, data = barnacle, kbin = 200,
#'                nboot = 100, smooth = "kernel")
#' 
#' # Nonparametric regression with interactions
#' fit2 <- frfast(DW ~ RC : F, data = barnacle, nboot = 100)
#' fit2
#' summary(fit2)
#' 
#' # using  splines
#' #fit2 <- frfast(DW ~ s(RC, by = F), data = barnacle,
#' #               nboot = 100, smooth = "splines", cluster = TRUE, ncores = 2)
#' #fit2
#' #summary(fit2)
#' 
#' 
#' # Allometric model
#' fit3 <- frfast(DW ~ RC, data = barnacle, model = "allo", nboot = 100)
#' summary(fit3)
#' 
#' # fit4 <- frfast(DW ~ RC : F, data = barnacle, model = "allo", nboot = 100)
#' # summary(fit4)
#' 
#' @useDynLib npregfast frfast_
#' @importFrom stats na.omit runif lm predict quantile
#' @importFrom mgcv interpret.gam gam predict.gam
#' @importFrom sfsmisc D1D2
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores
#' @importFrom foreach foreach %dopar%
#' 
#' @export




frfast <- function(formula, data = data, model = "np", smooth = "kernel", 
                   h0 = -1.0, h = -1.0, 
                   nh = 30, weights = NULL, kernel = "epanech", p = 3, 
                   kbin = 100, nboot = 500, rankl = NULL, ranku = NULL, 
                   seed = NULL, cluster = TRUE, ncores = NULL, ...){
  
  if(kernel == "gaussian")  kernel <- 3
  if (kernel == "epanech")   kernel <- 1
  if (kernel == "triang")    kernel <- 2
  
  if (missing(formula)) {
    stop("Argument \"formula\" is missing, with no default")
  }
  if (missing(data)) {
    stop("Argument \"data\" is missing, with no default")
  }
  if (!(kernel %in% 1:3)) {
    stop("Kernel not suported")
  }
  
  if (!(smooth %in% c("kernel", "splines"))) {
    stop("Smoother not suported")
  }
  
  #if(is.null(seed)) seed <- -1
  
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
  
  
  
  
  
  if (is.null(f)) f <- rep(1.0, n)
  etiquetas <- unique(f)
  nf <- length(etiquetas)
  
  if (model == "np") tmodel <- 1
  if (model == "allo") tmodel <- 2
  
  
  
  if (is.null(h0)) {
    h0 <- -1.0
  }
  if (is.null(h)) {
    h <- rep(-1.0, nf)
  }else{
    if (length(h) == 1) h <- rep(h, nf)
  }
  
  
  if (is.null(weights)) {
    weights <- rep(1.0, n)
  }else{
    if (sum(weights) <= 0 || any(weights < 0) || length(weights) != n) 
      stop("The specified weights are not correct")
  }  
  
  
  
  if (is.null(c2)) c2 <- matrix(as.double(-1.0), ncmax, nf) 
  if (is.null(rankl)) {
    rankl <- na.omit(as.vector(tapply(data[ ,varnames], f, min)))
  }else{
    if (length(rankl) == 1) rankl <- rep(rankl, nf)
  }
  if (is.null(ranku)) {
    ranku <- na.omit(as.vector(tapply(data[ ,varnames], f, max)))
  }else{
    if (length(ranku) == 1) ranku <- rep(ranku, nf)
  } 
  
  
  
  
  
  
  
  if (smooth != "splines") {
    
    ipredict2 <- 0
    
    umatrix <- matrix(runif(n*nboot), ncol = nboot, nrow = n)
    
    frfast  <- .Fortran("frfast_",
                        f = as.integer(f),
                        x = as.double(data[ ,varnames]),
                        y = as.double(data[ ,ffr$response]),
                        w = as.double(weights),
                        n = as.integer(n),
                        h0 = as.double(h0),
                        h = as.double(h),
                        c2 = as.integer(c2),
                        ncmax = as.integer(ncmax),
                        p = as.integer(p),
                        kbin = as.integer(kbin),
                        #fact = as.integer(c(1:nf)), 
                        fact = unique(as.integer(f)), 
                        nf = as.integer(nf),
                        nboot = as.integer(nboot),
                        xb = as.double(rep(-1.0, kbin)),
                        pb = array(rep(-1.0), c(kbin, 3, nf)),
                        li = array(as.double(-1.0), c(kbin, 3, nf)),
                        ls = array(as.double(-1.0), c(kbin, 3, nf)),
                        dif = array(as.double(-1.0), c(kbin, 3, nf, nf)),
                        difi = array(as.double(-1.0), c(kbin, 3, nf, nf)),
                        difs = array(as.double(-1.0), c(kbin, 3, nf, nf)),
                        tmodel = as.integer(tmodel), 
                        c = array(as.double(-1.0), c(3, nf)),
                        cs = array(as.double(-1.0), c(3, nf)),
                        ci = array(as.double(-1.0), c(3, nf)),
                        difc = array(as.double(-1.0), c(3, nf, nf)),
                        difcs = array(as.double(-1.0), c(3, nf, nf)),
                        difci = array(as.double(-1.0), c(3, nf, nf)),
                        pboot = array(as.double(-1.0), c(kbin, 3, nf, nboot)),
                        pcmin = as.double(rankl),
                        pcmax = as.double(ranku), 
                        cboot = array(as.double(-1.0), c(3, nf, nboot)), 
                        kernel = as.integer(kernel),
                        nh = as.integer(nh),
                        a = as.double(rep(-1, nf)),
                        ainf = as.double(rep(-1, nf)),
                        asup = as.double(rep(-1, nf)),
                        b = as.double(rep(-1, nf)),
                        binf = as.double(rep(-1, nf)),
                        bsup = as.double(rep(-1, nf)),
                        ipredict = as.integer(ipredict2),
                        predict = array(rep(-1.0), c(kbin, 3, nf)),
                        predictl = array(as.double(-1.0), c(kbin, 3, nf)),
                        predictu = array(as.double(-1.0), c(kbin, 3, nf)),
                        seed = as.integer(seed),
                        umatrix = as.double(umatrix)
    )
    
    if (tmodel != 2) {
      frfast$a <- NULL
      frfast$ainf <- NULL
      frfast$asup <- NULL
      frfast$b <- NULL
      frfast$binf <- NULL
      frfast$bsup <- NULL
      r2 <- NULL
    }
    
    #R-squared
    if (tmodel == 2) {
      yhat <- frfast$a * (frfast$x ^ frfast$b)
      rss <- sum( (frfast$y - yhat) ** 2 ) / (frfast$n - 2)
      tss <- sum(  (frfast$y - mean(frfast$y)) ** 2 ) / (frfast$n - 1)
      r2 <- 1 - (rss/tss)
    }
    
    frfast$pb[frfast$pb == -1] <- NA
    frfast$li[frfast$li == -1] <- NA
    frfast$ls[frfast$ls == -1] <- NA
    frfast$dif[frfast$dif == -1] <- NA
    frfast$difi[frfast$difi == -1] <- NA
    frfast$difs[frfast$difs == -1] <- NA
    frfast$c[frfast$c == -1] <- NA
    frfast$cs[frfast$cs == -1] <- NA
    frfast$ci[frfast$ci == -1] <- NA
    frfast$difc[frfast$difc == -1] <- NA
    frfast$difcs[frfast$difcs == -1] <- NA
    frfast$difci[frfast$difci == -1] <- NA
    frfast$pboot[frfast$pboot == -1] <- NA
    
    res <- list(x = frfast$xb,
                p = frfast$pb, 
                pl = frfast$li,
                pu = frfast$ls,
                diff = frfast$dif,
                diffl = frfast$difi,
                diffu = frfast$difs,
                nboot = frfast$nboot,
                n = frfast$n,
                dp = frfast$p,
                h = frfast$h,
                h0 = frfast$h0,
                fmod = frfast$f,
                xdata = as.vector(frfast$x),
                ydata = frfast$y,
                w = frfast$w,
                #fact=fact,  # Lo tuve que comentar pq me daba error
                kbin = frfast$kbin,
                nf = frfast$nf,
                max = frfast$c, #
                maxu = frfast$cs, 
                maxl = frfast$ci,
                diffmax = frfast$difc,
                diffmaxu = frfast$difcs,
                diffmaxl = frfast$difci,
                repboot = frfast$pboot,  
                rankl = frfast$pcmin,
                ranku = frfast$pcmax,
                nmodel = frfast$tmodel, 
                label = as.character(etiquetas),
                numlabel = unique(frfast$f),
                kernel = frfast$kernel, 
                a = frfast$a,
                al = frfast$ainf,
                au = frfast$asup,
                b = frfast$b,
                bl = frfast$binf,
                bu = frfast$bsup,
                name = c(ffr$response,varnames),
                formula = formula,
                nh = frfast$nh,
                r2 = r2,
                smooth = smooth,
                cluster = NA,
                call = match.call()
    )
    
  }else{
    
    mainfun <- function(formula, data, weights){
      
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
      
      # differences curves
      diff <- array(NA, dim = c(kbin, 3, nf, nf))
      if (nf > 1) {
        aux <- combn(1:nf, m = 2)
        for (j in 1:ncol(aux)) {
          diff[, , aux[1, j], aux[2, j]] <- p[, , aux[1, j]] - p[, , aux[2, j]]
        }    
      }
      
      
      
      # differences max
      diffmax <- array(NA, dim = c(3, nf, nf))
      if (nf > 1) {
        aux <- combn(1:nf, m = 2)
        for (j in 1:ncol(aux)) {
          diffmax[, aux[1, j], aux[2, j]] <- max[, aux[2, j]] - max[, aux[1, j]]
        }  
      }
      
      
      res <- list()
      res$p <- p
      res$diff <- diff
      res$max <- max
      res$diffmax <- diffmax  
      res$xgrid <- xgrid
      return(res)
    }
    
    
    
    
    
    
    
    res <- mainfun(formula, data = data, weights = weights)
    
    #  for(j in 1:nf){
    #   res$max[res$max[, j] == ranku[j], j] = NA
    #  }
    
    # bootstrap
    m <- gam(formula, weights = weights, data = data.frame(data, weights), ...)
    muhat <- as.vector(predict(m, type = "response"))
    err <- data[, ffr$response] - muhat
    err <- err - mean(err)
    yboot <- replicate(nboot, muhat + err *
                         sample(c(-sqrt(5) + 1, sqrt(5) + 1)/2, size = n,
                                replace = TRUE,
                                prob = c(sqrt(5) + 1, sqrt(5) - 1)/(2 * sqrt(5))))
    
    
    
    allboot <- foreach(i = 1:nboot) %dopar% {
      datab <- data
      datab$DW <- yboot[, i]
      aux <- mainfun(formula, data = data.frame(datab, weights), 
                     weights = weights, ...)
      return(aux)
    }
    
    
    
    pboot <- lapply(allboot, function(x){x$p})
    pboot <- array(unlist(pboot), dim = c(kbin, 3, nf, nboot))
    
    diffboot <- lapply(allboot, function(x){x$diff})
    diffboot <- array(unlist(diffboot), dim = c(kbin, 3, nf, nf, nboot))
    
    maxboot <- lapply(allboot, function(x){x$max})
    maxboot <- array(unlist(maxboot), dim = c(3, nf, nboot))
    
    diffmaxboot <- lapply(allboot, function(x){x$diffmax})
    diffmaxboot <- array(unlist(diffmaxboot), dim = c(3, nf, nf, nboot))
    
    
    # ci  p
    aux <- apply(pboot, c(3,2,1), 
                 function(x){quantile(x, probs = c(0.025), na.rm = TRUE)})
    aux2 <- apply(pboot, c(3,2,1), 
                  function(x){quantile(x, probs = c(0.975), na.rm = TRUE)})
    
    pl <- array(NA, dim = c(kbin, 3, nf))
    pu <- array(NA, dim = c(kbin, 3, nf))
    for (i in  1:nf) {
      pl[, , i] <- t(aux[i, , ])
      pu[, , i] <- t(aux2[i, , ])
    }
    
    
    # ci  diff 
    diffl <- array(NA, dim = c(kbin, 3, nf, nf))
    diffu <- array(NA, dim = c(kbin, 3, nf, nf))
    if (nf > 1) {
      com <- combn(1:nf, m = 2)
      for (j in 1:ncol(com)) {
        diffl[, , com[1,j], com[2,j] ] <- t(apply(diffboot[, , com[1,j], com[2,j], ],c(2,1), 
                                                  function(x){quantile(x, probs = c(0.025), na.rm = TRUE)}))
      }
      for (j in 1:ncol(com)) {
        diffu[, , com[1,j], com[2,j] ] <- t(apply(diffboot[, , com[1,j], com[2,j], ],c(2,1),
                                                  function(x){quantile(x, probs = c(0.975), na.rm = TRUE)}))
      }
    }
    
    
    
    # ci  max
    aux <- apply(maxboot, c(2,1), 
                 function(x){quantile(x, probs = c(0.025), na.rm = TRUE)})
    aux2 <- apply(maxboot, c(2,1), 
                  function(x){quantile(x, probs = c(0.975), na.rm = TRUE)})
    
    maxl <- t(aux)
    maxu <- t(aux2)
    
    for(j in 1:nf){
      maxl[maxl[, j] == ranku[j], j] = NA
      maxu[maxu[, j] == ranku[j], j] = NA
    }
    
    # ci  diffmax FALTA
    diffmaxl <- array(NA, dim = c(3, nf, nf))
    diffmaxu <- array(NA, dim = c(3, nf, nf))
    if (nf > 1) {
      com <- combn(1:nf, m = 2)
      for (j in 1:ncol(com)) {
        diffmaxl[, com[1,j], com[2,j] ] <- t(apply(diffmaxboot[, com[1,j], com[2,j], ], c(1), 
                                                   function(x){quantile(x, probs = c(0.025), na.rm = TRUE)}))
        diffmaxu[, com[1,j], com[2,j] ] <- t(apply(diffmaxboot[, com[1,j], com[2,j], ], c(1), 
                                                   function(x){quantile(x, probs = c(0.975), na.rm = TRUE)}))
      }
    }
    #  maxl <- t(aux)
    #  maxu <- t(aux2)
    #diffmaxu <- NA
    #diffmaxl <- NA
    #}
    
    
    res <- list(x = res$xgrid,
                p = res$p, 
                pl = pl,
                pu = pu,
                diff = res$diff,
                diffl = diffl,
                diffu = diffu,
                nboot = nboot,
                n = n,
                dp = NA,
                h = NA,
                h0 = NA,
                fmod = f,
                xdata = data[, varnames],
                ydata = data[, ffr$response],
                w = weights,
                #fact=fact,  # Lo tuve que comentar pq me daba error
                kbin = kbin,
                nf = nf,
                max = res$max, #
                maxu = maxu, 
                maxl = maxl,
                diffmax = res$diffmax,
                diffmaxu = diffmaxu,
                diffmaxl = diffmaxl,
                repboot = pboot,  
                rankl = rankl,
                ranku = ranku,
                nmodel = tmodel, 
                label = as.character(etiquetas),
                numlabel = unique(f),
                kernel = kernel, 
                a = NA,
                al = NA,
                au = NA,
                b = NA,
                bl = NA,
                bu = NA,
                name = c(ffr$response,varnames),
                formula = formula,
                nh = nh,
                r2 = NA,
                smooth = smooth,
                cluster = cluster,
                ncores = ncores,
                call = match.call()
    )
  }
  
  
  class(res) <- "frfast"
  return(res)
}


