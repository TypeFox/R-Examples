###### mean function ########################################################
skewhypMean <- function (mu = 0, delta = 1, beta = 1, nu = 1,
                  param = c(mu,delta,beta,nu)) {

    #check parameters
    parResult <- skewhypCheckPars(param)
    case <- parResult$case
    errMessage <- parResult$errMessage
    if(case == "error") stop(errMessage)
    mu <- param[1]
    delta <- param[2]
    beta <- param[3]
    nu <- param[4]

    skewhypMean <- mu + (beta*delta^2)/(nu-2)

    return(skewhypMean)

}

###### variance function #####################################################
skewhypVar<- function (mu = 0, delta = 1, beta = 1, nu = 1,
                     param = c(mu,delta,beta,nu)) {

    #check parameters
    parResult <- skewhypCheckPars(param)
    case <- parResult$case
    errMessage <- parResult$errMessage
    if(case == "error") stop(errMessage)
    mu <- param[1]
    delta <- param[2]
    beta <- param[3]
    nu <- param[4]

    if (nu <= 4)stop ("variance does not exist when nu <= 4")

    skewhypVar<-((2*beta^2*delta^4)/((nu-2)^2*(nu-4)))+(delta^2)/(nu-2)

    return(skewhypVar)
}
####### skewness function ##################################################
skewhypSkew<-function(mu = 0, delta = 1, beta = 1, nu = 1,
                   param = c(mu,delta,beta,nu)){

    #check parameters
    parResult <- skewhypCheckPars(param)
    case <- parResult$case
    errMessage <- parResult$errMessage
    if(case == "error") stop(errMessage)
    mu <- param[1]
    delta <- param[2]
    beta <- param[3]
    nu <- param[4]

    if(nu <= 6)stop("skewness does not exist when nu <= 6")

    skewhypSkew<- ((2*(nu-4)^0.5*beta*delta)/
                   ((2*beta^2*delta^2+(nu-2)*(nu-4))^(3/2)))*
                       (3*(nu-2)+(8*beta^2*delta^2)/(nu-6))

    return(skewhypSkew)

}

###### kurtosis function ####################################################
skewhypKurt <- function (mu = 0, delta = 1, beta = 1, nu = 1,
                         param = c(mu,delta,beta,nu)) {

    #check parameters
    parResult <- skewhypCheckPars(param)
    case <- parResult$case
    errMessage <- parResult$errMessage
    if (case == "error") stop(errMessage)
    mu <- param[1]
    delta <- param[2]
    beta <- param[3]
    nu <- param[4]

    if (nu <= 8) stop("kurtosis does not exist when nu <= 8")

    skewhypKurt <- (6/((2*beta^2*delta^2+(nu-2)*(nu-4))^2))*
        ((nu-2)^2*(nu-4)+(16*beta^2*delta^2*(nu-2)*(nu-4))/(nu-6)
         +(8*beta^4*delta^4*(5*nu-22))/(nu-6)*(nu-8))

    return(skewhypKurt)
}

###### mode function ########################################################
skewhypMode <- function(mu = 0, delta = 1, beta = 1, nu = 1,
                  param = c(mu,delta,beta,nu),
                  tolerance = .Machine$double.eps ^ 0.5){


    #check parameters
    parResult <- skewhypCheckPars(param)
    case <- parResult$case
    errMessage <- parResult$errMessage
    if(case == "error") stop(errMessage)
    mu <- param[1]
    delta <- param[2]
    beta <- param[3]
    nu <- param[4]

    if (abs(beta) > tolerance) { #beta != 0

        modeFun <- function(x){
            dskewhyp(x,param=param,log=TRUE)
        }
        start <- mu
        options(warn=-1)
        opt <- optim(start, modeFun, control = list(fnscale = -1,
                                           maxit = 1000, method = "BFGS"))
        ifelse (opt$convergence == 0, distMode <- opt$par, distMode <- NA)
    }
    else{ #beta = 0
        distMode <- mu
    }

    return(distMode)
}



