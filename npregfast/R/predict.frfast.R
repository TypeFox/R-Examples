#' Prediction from fitted \code{frfast} model
#' @description Takes a fitted \code{frfast} object  and produces predictions 
#' (with their 95\% confidence intervals) from a fitted model with 
#' interactions or without interactions.
#' @param object A fitted \code{frfast} object as produced by \code{frfast()}.
#' @param newdata A data frame containing the values of the model covariates
#' at which predictions are required.  If newdata is provided, then it should 
#' contain all the variables needed for prediction: a warning is 
#' generated if not.
#' @param fac Factor's level to take into account. By default is \code{NULL}.
#' @param der Number which determines any inference process. By default 
#' \code{der} is \code{NULL}. If this term is \code{0}, 
#' the function returns the initial estimate. If it is \code{1} or \code{2}, 
#' it is designed for the first or second derivative, respectively.
#' @param seed Seed to be used in the bootstrap procedure.
#' @param \ldots Other options.
#' 
#' @return \code{predict.frfast} computes and returns a list containing 
#' predictions of the estimates, first and second derivative, 
#' with their 95\% confidence intervals.
#' 
#'@author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#'@examples
#' library(npregfast)
#' data(barnacle)
#' 
#' # Nonparametric regression without interactions
#' fit <- frfast(DW ~ RC, data = barnacle, nboot = 100)
#' nd <- data.frame(RC = c(10, 14, 18))
#' predict(fit, newdata = nd)
#' 
#' # Nonparametric regression with interactions
#' # fit2 <- frfast(DW ~ RC : F, data = barnacle, nboot = 100)
#' # nd2 <- data.frame(RC = c(10, 15, 20))
#' # predict(fit2, newdata = nd2)
#' # predict(fit2, newdata = nd2, der = 0, fac = "barca")
#'  
#' @useDynLib npregfast frfast_
#' @importFrom stats na.omit runif
#' @export

predict.frfast <- function(object = model, newdata, fac = NULL, der = NULL, 
                           seed = NULL, ...) {
  
  model <- object
  
  if (missing(newdata)) {
    stop("Argument \"newdata\" is missing, with no default")
  }
  
  if (length(fac) > 1) {
    stop("Argument \"fac\" have to be a length-one vector")
  }
  
  if (length(der) > 1) {
    stop("Argument \"der\" have to be a length-one vector")
  }
  
  
  if (!is.null(fac) & model$nf == 1) {
    stop("Argument \"fac\" not supported. 
         There is not factor in the model.")
  }
  
  if (!is.null(fac) & !isTRUE(fac %in% model$label)) {
    stop("\"",paste(fac),"\" is not a factor's level.")
  }
  
  if (!is.null(der) & !isTRUE(der %in% c(0, 1, 2))) {
    stop("",paste(der)," is not a r-th derivative implemented, only 
         permitted 0, 1 or 2.")
  }
  
 # if(is.null(seed)) seed <- -1
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (isTRUE(model$cluster) & model$smooth == "splines") {
      num_cores <- model$ncores
    registerDoParallel(cores = num_cores)
  }
  
  
  if (model$smooth != "splines") {
  
  newdata <- newdata[, 1]
  len <- length(newdata)
  f <- model$fmod
  f <- c(f, rep(1, length(newdata)))  
  x <- model$xdata
  x <- c(x, newdata)
  y <- model$ydata
  y <- c(y, rep(-1.0, length(newdata)))  
  n <- length(x)
  h0 <- model$h0
  h <- model$h
  weights <- c(model$w, rep(0, length(newdata)))
  p <- model$dp
  kbin <- model$kbin
  ncmax <- 5
  ikernel <- 1
  iopt <- 1
  nboot <- model$nboot
  rankl <- model$rankl
  ranku <- model$ranku
  kernel <- model$kernel
  nh <- model$nh
  nf <- model$nf
  c2 <- matrix(as.double(-1.0), ncmax, nf) 
  tmodel <- model$nmodel
  ipredict2 <- 1
  
  umatrix <- matrix(runif(n*nboot), ncol = nboot, nrow = n)
  
  frfast  <- .Fortran("frfast_",
                      f      = as.integer(f),
                      x      = as.double(x),
                      y      = as.double(y),
                      w      = as.double(weights),
                      n      = as.integer(n),
                      h0      = as.double(h0),
                      h      = as.double(h),
                      c2     = as.integer(c2),
                      ncmax  = as.integer(ncmax),
                      p      = as.integer(p),
                      kbin   = as.integer(kbin),
                      fact   = as.integer(c(1:nf)),
                      nf     = as.integer(nf),
                      nboot  = as.integer(nboot),
                      xb     = as.double(rep(-1.0, kbin)),
                      pb     = array(rep(-1.0), c(kbin, 3, nf)),
                      li     = array(as.double(-1.0), c(kbin, 3, nf)),
                      ls     = array(as.double(-1.0), c(kbin, 3, nf)),
                      dif    = array(as.double(-1.0), c(kbin, 3, nf, nf)),
                      difi   = array(as.double(-1.0), c(kbin, 3, nf, nf)),
                      difs	  = array(as.double(-1.0), c(kbin, 3, nf, nf)),
                      tmodel  = as.integer(tmodel), 
                      c 	  = array(as.double(-1.0), c(3, nf)),
                      cs	  = array(as.double(-1.0), c(3, nf)),
                      ci	  = array(as.double(-1.0), c(3, nf)),
                      difc	  = array(as.double(-1.0), c(3, nf, nf)),
                      difcs  = array(as.double(-1.0), c(3, nf, nf)),
                      difci  = array(as.double(-1.0), c(3, nf, nf)),
                      pboot = array(as.double(-1.0), c(kbin, 3, nf, nboot)),
                      pcmin = as.double(rankl),
                      pcmax = as.double(ranku), 
                      cboot = array(as.double(-1.0), c(3, nf, nboot)), 
                      kernel = as.integer(kernel),
                      nh    = as.integer(nh),
                      a     = as.double(rep(-1.0, nf)),
                      ainf  = as.double(rep(-1.0, nf)),
                      asup  = as.double(rep(-1.0, nf)),
                      b     = as.double(rep(-1.0, nf)),
                      binf  = as.double(rep(-1.0, nf)),
                      bsup  = as.double(rep(-1.0, nf)),
                      ipredict = as.integer(ipredict2),
                      predict = array(as.double(-1.0), c(n, 3, nf)),
                      predictl = array(as.double(-1.0), c(n, 3, nf)),
                      predictu = array(as.double(-1.0), c(n, 3, nf)),
                      seed = as.integer(seed),
                      umatrix = as.double(umatrix)
  )
  
  
  frfast$predict[frfast$predict == -1] <- NA
  frfast$predictl[frfast$predictl == -1] <- NA
  frfast$predictu[frfast$predictu == -1] <- NA
  
  }else{
    
    
    formula <- model$formula
    ffr <- interpret.gam(formula)
    varnames <- ffr$pred.names[1]
    namef <- ffr$pred.names[2]
    response <- ffr$response
    data <- data.frame(model$xdata, model$ydata, model$fmod)
    names(data) <- c(paste(varnames), paste(response), paste(namef))
    weights <- model$w
    newdata <- newdata[, 1]
    len <- length(newdata)
  
    
    mainfun <- function(formula, data, weights){
      # grid
      xgrid <- seq(min(data[ ,varnames]), max(data[ ,varnames]), length.out = model$kbin)
      xgrid <- c(xgrid, newdata)
      newd <- expand.grid(xgrid, unique(model$fmod))
      names(newd) <- ffr$pred.names
      
      # estimations
      p <- array(NA, dim = c(model$kbin + length(newdata), 3, model$nf))
      m <- gam(formula, weights = weights, data = data.frame(data, weights), ...)
      muhat <- as.vector(predict(m, newdata = newd, type = "response"))
      p[, 1, 1:model$nf] <- muhat
      #d1 <- apply(p, 3, function(z){D1ss(x = xgrid, y = z[, 1])})
      d1 <- apply(p, 3, function(z){D1D2(x = xgrid, y = z[, 1], deriv = 1)$D1})
      p[, 2, 1:model$nf] <- as.vector(d1)
      #d2 <- apply(p, 3, function(z){D2ss(x = xgrid, y = z[, 1])$y})
      d2 <- apply(p, 3, function(z){D1D2(x = xgrid, y = z[, 1], deriv = 2)$D2})
      p[, 3, 1:model$nf] <- as.vector(d2)
      return(p)
    }
      
    
    
    
    
    res <- mainfun(formula, data = data, weights = weights)
    p <- res

    
    # bootstrap
    m <- gam(formula, weights = weights, data = data.frame(data, weights), ...)
    muhat <- as.vector(predict(m, type = "response"))
    err <- data[, ffr$response] - muhat
    err <- err - mean(err)
    yboot <- replicate(model$nboot, muhat + err *
                         sample(c(-sqrt(5) + 1, sqrt(5) + 1)/2, size = model$n,
                                replace = TRUE,
                                prob = c(sqrt(5) + 1, sqrt(5) - 1)/(2 * sqrt(5))))
    
    
    
    allboot <- foreach(i = 1:model$nboot) %dopar% {
      datab <- data
      datab$DW <- yboot[, i]
      aux <- mainfun(formula, data = data.frame(datab, weights), 
                     weights = weights, ...)
      return(aux)
    }
    
    #pboot <- lapply(allboot, function(x){x$p})
    pboot <- array(unlist(allboot), dim = c(model$kbin + length(newdata), 3, model$nf, model$nboot))
    
    
    
    # ci  p
    aux <- apply(pboot, c(3,2,1), 
                 function(x){quantile(x, probs = c(0.025), na.rm = TRUE)})
    aux2 <- apply(pboot, c(3,2,1), 
                  function(x){quantile(x, probs = c(0.975), na.rm = TRUE)})
    
    pl <- array(NA, dim = c(model$kbin + length(newdata), 3, model$nf))
    pu <- array(NA, dim = c(model$kbin + length(newdata), 3, model$nf))
    for (i in  1:model$nf) {
      pl[, , i] <- t(aux[i, , ])
      pu[, , i] <- t(aux2[i, , ])
    }
    
    
  } # close else
  
  
  
  
  if (model$smooth != "splines") {
   ii <- c(rep(FALSE, n - length(newdata)), rep(TRUE, length(newdata)))
   aux <- frfast$predict
   auxl <- frfast$predictl
   auxu <- frfast$predictu
  }else{
    ii <- c(rep(FALSE, model$kbin) , rep(TRUE, length(newdata)))
    aux <- p
    auxl <- pl
    auxu <- pu
  }
  
nf <- model$nf
  
  if (nf == 1) {
    for (k in 1:nf) {
      if (is.null(der)) {
        der <- c(0, 1, 2)
      }
      cont <- 0
      der <- der + 1
      res <- array(data = NA, dim = c(len, 3, length(der)))
      for (j in der) {
        cont <- cont + 1
        res[, 1, cont] <- aux[, , 1][ii, j]
        res[, 2, cont] <- auxl[, , 1][ii, j]
        res[, 3, cont] <- auxu[, , 1][ii, j]
      }
    }
    colnames(res) <- c("Pred", "Lwr", "Upr")
    if (length(der) == 3) {
      res <- list(Estimation = res[, , 1], First_deriv = res[, , 2], 
                  Second_deriv = res[, , 3])
    } else {
      
      if (der == 1) {
        res <- list(Estimation = res[, , 1])
      }
      if (der == 2) {
        res <- list(First_deriv = res[, , 1])
      }
      if (der == 3) {
        res <- list(Second_deriv = res[, , 1])
      }
    }
    
    class(res) <- "predict.frfast"
    return(res)
  } else {
    if (is.null(fac)) {
      fac <- unique(model$fmod)
    }else{
      fac <- which(model$label == fac)
    }
    factores <- c()
    resul <- vector("list", length = length(fac))
    zz <- 1
    for (k in 1:length(fac)) {
     if (model$smooth != "splines") {
       factores[zz] <- paste("Level_", model$label[frfast$fact[k]], sep = "")
     }else{
       factores[zz] <- paste("Level_", model$label[k], sep = "")
     }
      zz <- zz + 1
    }
    names(resul) <- factores
    if (is.null(der)) {
      der <- c(0, 1, 2)
    }
    der <- der + 1
    zz <- 1
    for (k in 1:length(fac)) {
      cont <- 0
      res <- array(data = NA, dim = c(len, 3, length(der)))
      for (j in der) {
        cont <- cont + 1
        res[, 1, cont] <- aux[, , k][ii, j]
        res[, 2, cont] <- auxl[, , k][ii, j]
        res[, 3, cont] <- auxu[, , k][ii, j]
      }
      colnames(res) <- c("Pred", "Lwr", "Upr")
      # res<-list(Estimation=res[,,1], First_deriv=res[,,2], Second_deriv=res[,,3])
      if (length(der) == 3) {
        res <- list(Estimation = res[, , 1], First_deriv = res[, , 2], 
                    Second_deriv = res[, , 3])
      } else {
        
        if (der == 1) {
          res <- list(Estimation = res[, , 1])
        }
        if (der == 2) {
          res <- list(First_deriv = res[, , 1])
        }
        if (der == 3) {
          res <- list(Second_deriv = res[, , 1])
        }
      }
      
      resul[[zz]] <- res
      zz <- zz + 1
    }
  }
    class(resul) <- "predict.frfast"
    return(resul)
  }
  


