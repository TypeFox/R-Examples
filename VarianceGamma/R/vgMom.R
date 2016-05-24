## calculate VG moments (mu moments, central moments, raw moments or others)
## if only momType is defined, then calculate moments according to momType.
## if only about is defined, then calculate moment according to about.
## if users define both momType and about(contradictive or not), then
## the about value would always overwrites momType. Thus, calculate momoents
## about "about".

vgMom <- function(order, vgC = 0, sigma = 1, theta = 0, nu = 1,
  param = c(vgC,sigma,theta,nu), momType = "raw", about = 0) {

    ## check order is whole number
    if (!is.wholenumber(order)){
      stop("Order must be a whole number")
    }
    if ((order < 0)) {
      stop("Order must be positive")
    }

    ## check momType
    momType <- as.character(momType)
    momType <- tolower(momType)
    if (momType != "raw" & momType != "central" & momType != "mu") {
      stop ("Unrecognised moment type")
    }

    ## check parameters
    parResult <- vgCheckPars(param = param)
    case <- parResult$case
    errMessage <- parResult$errMessage
    if (case == "error"){
      stop(errMessage)
    } else if (case =="normal"){
        if (order == 0) {
          mom <- 1
        } else {
            vgC <- param[1]
            sigma <- param[2]
            theta <- param[3]
            nu <- param[4]
            ## Transform vg parameters to GH parameters
            ghParam <- vgChangePars(1, 4, param = param)
            lambda <- ghParam[1]
            alpha <- ghParam[2]
            beta <- ghParam[3]
            ## calculate mu moments
            muMom <- rep (NA,order)
            for (i in 1:order) {
              a <- momRecursion(order = i)
              coeff <- a$a
              betaPow <- a$M
              gammaPow <- a$L
              lengthGammaPow <- length(gammaPow)
              muM <- coeff*beta^betaPow*
                    (sapply(gammaPow[1]:gammaPow[lengthGammaPow], gammaRawMom,
                                  shape = lambda, rate = (alpha^2 - beta^2)/2))
              muMom[i] <- sum(muM)
            }

            if (about != 0) {
              mom <- momChangeAbout(order = order, oldMom = muMom,
                                    oldAbout = vgC, newAbout = about)
            } else {
                if (momType == "mu") {
                  mom = muMom[order]
                } else if (momType == "raw") {
                  about = 0
                  mom <- momChangeAbout(order = order, oldMom = muMom,
                                        oldAbout = vgC, newAbout = about)
                } else if (momType == "central") {
                  about = vgMean (param = param)
                  mom <- momChangeAbout(order = order, oldMom = muMom,
                                        oldAbout = vgC, newAbout = about)
                }
            }

        }

    }
    return(mom)
}
