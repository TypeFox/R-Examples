etas_normal0 <-
function(data, evalpts, params, fixedparams, TT = NA, tplus = FALSE)
{
#   spatial effect like a normal density function but with
#   only exponential term
#   major and minor axes in direction of x and y but with
#   different scaling factors
#   message("NOTE: 'etas.normal0' and its methods are still under development and may change.")
    if (is.null(data)) data <- cbind(time = Inf, magnitude = 0)
    #------------------------------------------------------------
    if(any(is.na(TT))) {
        lambda <- function(oneevalpt, data, params, tplus)
        {
            #  Note: lambda(t) only evaluates one point (t scalar)
            #-----------------------------------
            #   Parameters
            mu <- params[1]
            A <- params[2]
            alpha <- params[3]
            CC <- params[4]
            P <- params[5]
            dx <- params[6]
            dy <- params[7]
            beta <- params[8]
            #-----------------------------------
            tj <- oneevalpt["time"]
            ti <- data[, "time"]
            if (!tplus) use <- (ti < tj)
            else use <- (ti <= tj)
            if(any(use)){
                #-----------------------------------
                #   History
                Mi <- data[use, "magnitude"]
                ti <- data[use, "time"]
                xi <- data[use, "longitude"]
                yi <- data[use, "latitude"]
                #-----------------------------------
                #   Evaluation Point
                xj <- oneevalpt["longitude"]

                yj <- oneevalpt["latitude"]
                magns <- exp(alpha*Mi)
                omori <- (1+(tj-ti)/CC)^(-P)
                spatial <- exp( -(xj-xi)^2/(2*dx*magns) - (yj-yi)^2/(2*dy*magns) )
                ci <- mu + A*sum(exp(beta*Mi)*omori*spatial)
            } else ci <- mu
            return(ci)
        }
        #---------------------------------------
#       if (exists("SNOWcluster")){
#           if (inherits(SNOWcluster, "cluster"))
#               ci <- parRapply(SNOWcluster, evalpts, lambda, data=data,
#                               params=params, tplus=tplus)
#           else stop("ERROR: Object SNOWcluster has wrong class")
#       } else ci <- apply(evalpts, 1, lambda, data=data,
#                          params=params, tplus=tplus)
        ci <- apply(evalpts, 1, lambda, data=data,
                    params=params, tplus=tplus)
    }
    else {
        #  if TT is not missing, then just evaluate the integral term
        mu <- params[1]
        A <- params[2]
        alpha <- params[3]
        CC <- params[4]
        P <- params[5]
        dx <- params[6]
        dy <- params[7]
        beta <- params[8]
        X1 <- fixedparams[1]
        X2 <- fixedparams[2]
        Y1 <- fixedparams[3]
        Y2 <- fixedparams[4]
        area <- (X2-X1)*(Y2-Y1)
        times <- data[, "time"]
        S <- c(NA, NA)
        for(k in 1:2) {
            use <- (times < TT[k])
            if(any(use)) {
                magns <- exp((alpha)*data[use,"magnitude"])
                sqrtmagns <- sqrt(magns)
                Fx <- pnorm(X2-data[use,"longitude"], sd=sqrt(dx)*sqrtmagns) -
                      pnorm(X1-data[use,"longitude"], sd=sqrt(dx)*sqrtmagns)
                Fy <- pnorm(Y2-data[use,"latitude"], sd=sqrt(dy)*sqrtmagns) -
                      pnorm(Y1-data[use,"latitude"], sd=sqrt(dy)*sqrtmagns)
                spatialint <- 2*pi*magns*sqrt(dx*dy)*Fx*Fy
                if(P!=1) omoriint <- (1 - (1 + (TT[k] -
                                 times[use])/CC)^(-P + 1))/(P - 1)
                else omoriint <- log(1 + (TT[k] - times[use])/CC)
                S[k] <- sum(exp(beta*data[use,"magnitude"])*spatialint*omoriint)
            }
            else S[k] <- 0
        }
        ci <- mu*(TT[2]-TT[1])*area + CC*A*(S[2]-S[1])
    }
    names(ci) <- NULL
    return(ci)
}


