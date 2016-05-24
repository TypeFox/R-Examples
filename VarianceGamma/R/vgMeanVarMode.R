vgMean <- function (vgC = 0, sigma = 1, theta = 0, nu = 1,
                    param = c(vgC,sigma,theta,nu)) {

  #check parameters
  parResult <- vgCheckPars(param = param)
  case <- parResult$case
  errMessage <- parResult$errMessage
  if (case == "error"){
    stop(errMessage)
  }
  vgC <- param[1]
  sigma <- param[2]
  theta <- param[3]
  nu <- param[4]

  distMean <- vgC + theta
  return(distMean)
}

vgVar <- function (vgC = 0, sigma = 1, theta = 0, nu = 1,
                   param = c(vgC,sigma,theta,nu)) {

  #check parameters
  parResult <- vgCheckPars(param = param)
  case <- parResult$case
  errMessage <- parResult$errMessage
  if (case == "error"){
    stop(errMessage)
  }
  vgC <- param[1]
  sigma <- param[2]
  theta <- param[3]
  nu <- param[4]

  distVar <- sigma^2 + theta^2*nu
  return(distVar)
}

vgSkew <- function (vgC = 0, sigma = 1, theta = 0, nu = 1,
                    param = c(vgC,sigma,theta,nu)) {

  #check parameters
  parResult <- vgCheckPars(param = param)
  case <- parResult$case
  errMessage <- parResult$errMessage
  if (case == "error"){
    stop(errMessage)
  }
  vgC <- param[1]
  sigma <- param[2]
  theta <- param[3]
  nu <- param[4]

  distSkew <- (2*theta^3*nu^2 + 3*sigma^2*theta*nu)/
                 ((theta^2*nu + sigma^2)^(3/2))
  return(distSkew)
}

vgKurt <- function (vgC = 0, sigma = 1, theta = 0, nu = 1,
                    param = c(vgC,sigma,theta,nu)) {

  #check parameters
  parResult <- vgCheckPars(param = param)
  case <- parResult$case
  errMessage <- parResult$errMessage
  if (case == "error"){
    stop(errMessage)
  }
  vgC <- param[1]
  sigma <- param[2]
  theta <- param[3]
  nu <- param[4]

  distKurt <- 3 + (3*sigma^4*nu + 12*sigma^2*theta^2*nu^2 + 6*theta^4*nu^3)/
                (theta^2*nu+sigma^2)^2
  return(distKurt)
}

vgMode <- function (vgC = 0, sigma = 1, theta = 0, nu = 1,
                    param = c(vgC,sigma,theta,nu)) {

    ##check parameters
    parResult <- vgCheckPars(param = param)
    case <- parResult$case
    errMessage <- parResult$errMessage
    if (case == "error"){
        stop(errMessage)
    }
    vgC <- param[1]
    sigma <- param[2]
    theta <- param[3]
    nu <- param[4]
    if (nu >= 2){
        distMode <- vgC
    } else {
        modeFun <- function(x) {
            log(dvg(x, param = param, log = FALSE))
        }
        start <- vgC
        optResult <- optim(start, modeFun,
                           control = list(fnscale = -1, maxit = 1000),
                           method = "BFGS")
        if (optResult$convergence == 0) {
            distMode <- optResult$par
        } else {
            distMode <- NA
        }
    }
    return(distMode)
}

