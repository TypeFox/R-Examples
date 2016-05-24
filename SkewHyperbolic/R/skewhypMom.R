############Calculate Skew Hyperbolic Moments############################

skewhypMom <- function(order, mu = 0, delta = 1, beta = 1, nu = 1,
                       param = c(mu, delta, beta, nu), momType = "raw",
                       about = 0){

    #check order
    if(!is.wholenumber(order)) stop("Order must be a whole number")
    if( order < 0) stop("Order must be positive")

    #check momType
    momType <- as.character(momType)
    momType <- tolower(momType)
    if(momType != "raw" & momType != "central" & momType != "mu")
        stop("Unrecognised moment type")

    #check parameters
    parResult <- skewhypCheckPars(param)
    case <- parResult$case
    errMessage <- parResult$errMessage

    if(case == "error") stop(errMessage)

    if(order == 0){
        mom <- 1
    }else{
        mu <- param[1]
        delta <- param[2]
        beta <- param[3]
        nu <- param[4]

        #calculate the mu moments
        muMom <- rep(NA, order)
        for(i in 1:order){
            a <- momRecursion(order = i)
            coeff <- a$a
            betaPow <- a$M
            gammaPow <- a$L
            lengthGammaPow <- length(gammaPow)
            muM <- coeff*beta^betaPow*
                (sapply(-gammaPow[1]:-gammaPow[lengthGammaPow],
                        gammaRawMom,
                        shape = (nu/2),
                        scale = (2/(delta^2))  ))#should this be scale or rate
            muMom[i] <- sum(muM)
        }
        if(about != 0){
            mom <- momChangeAbout(order = order, oldMom = muMom,
                                  oldAbout = mu, newAbout = about)
        }else{
            if(momType == "mu"){
                mom = muMom[order]
            }else if(momType == "raw"){
                about = 0
                mom <- momChangeAbout(order = order, oldMom = muMom,
                                      oldAbout = mu, newAbout = about)
            }else if(momType == "central"){
                about = skewhypMean(param = param)
                mom <- momChangeAbout(order = order, oldMom = muMom,
                                      oldAbout = mu, newAbout = about)
            }
        }
        return(mom)
    }
}

####### TO DO ###################################

#if not defined return NA
#momRecursion function - check the details of this....

