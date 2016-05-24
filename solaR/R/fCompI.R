fCompI <- function(sol, compD, G0I,
                   corr = 'none', f,
                   filterG0 = TRUE){
  
    ##Time indexes
    if (class(sol)=='Sol') {
        solI <- sol@solI
        mtch <- sol@match
        sample <- sol@sample
    } else { ##sol is a zoo, for example, produced by fSolI
        solI <- sol
        mtch <- attr(sol, 'match')
        sample <- attr(sol, 'sample')
    }
    indSol <- index(solI)
    ## Retrieve some variables from solI
    rd <- coredata(solI$rd)
    rg <- coredata(solI$rg)
    aman <- coredata(solI$aman)
    Bo0 <- coredata(solI$Bo0)

    ## If instantaneous values are not provided, compD is used instead.
    if (missing(G0I)) { 

        comp.rep <- data.frame(compD)[mtch, c('Ktd', 'G0d', 'D0d', 'B0d')]    

        ## Daily irradiation components
        D0d <- comp.rep$D0d
        G0d <- comp.rep$G0d
        B0d <- comp.rep$B0d

        ## Daily profile using Liu and Jordan, Collares-Pereira and Rabl
        ## proposals.
        D0 <- D0d * rd
        G0 <- G0d * rg
        ## This method may produce diffuse irradiance higher than global
        ## irradiance.
        G0 <- pmax(G0, D0, na.rm=TRUE)
        B0 <- G0 - D0

        ## Negative values are set to NA
        neg <- (B0 < 0)| (D0 < 0) | (G0 <0)
        is.na(G0) <- neg
        is.na(D0) <- neg
        is.na(B0) <- neg
    
        ## Daily profiles are scaled to keep daily irradiation values
        day <- truncDay(indSol)

        D0dCP <- ave(D0, day, FUN=function(x) P2E(x, sample))
        G0dCP <- ave(G0, day, FUN=function(x) P2E(x, sample))
        B0dCP <- ave(B0, day, FUN=function(x) P2E(x, sample))

        D0 <- D0 * D0d/D0dCP
        G0 <- G0 * G0d/G0dCP
        B0 <- B0 * B0d/B0dCP

        kt <- G0/Bo0
        fd <- D0/G0
  
    } else { ## Use instantaneous values if provided through G0I
    
        if (corr!='none'){ 
            if (class(G0I) == 'Meteo') {
                G0 <- coredata(getG0(G0I))
            } else {                       ## G0I is a zoo
                if (NCOL(G0I)>1) {         ## multivariable
                    G0 <- coredata(G0I$G0) # Only G0 is needed
                } else {
                    G0 <- coredata(G0I)
                }
            }                                 
            ## Filter values: surface irradiation must be lower than
            ## extraterrestial; 
            if (isTRUE(filterG0)) is.na(G0) <- (G0 > Bo0)

            kt <- G0/Bo0
    
            ## Fd-Kt correlation
            fd <- switch(corr,
                         EKDh = FdKtEKDh(kt),
                         CLIMEDh = FdKtCLIMEDh(kt),
                         BRL = FdKtBRL(kt, sol), 
                         user = f(kt, sol),      
                         stop('Wrong descriptor of the correlation fd-kt.'))
            D0 <- fd * G0
            B0 <- G0 - D0

        } else { ##corr=='none': G0 is a zoo with G0, D0 and B0

            if (class(G0I) == 'Meteo') {
                IrrData <- getData(G0I)
            } else {                    #G0I is a zoo
                IrrData <- G0I
            }
            D0 <- coredata(IrrData$D0)
            B0 <- coredata(IrrData$B0)
            G0 <- coredata(IrrData$G0)
            ## Filter values: surface irradiation must be lower than
            ## extraterrestial; 
            if (isTRUE(filterG0)) is.na(G0) <- is.na(D0) <- is.na(B0) <- (G0 > Bo0)
      
            kt <- G0/Bo0
            fd <- D0/G0
        }
    }
    ## Values outside sunrise-sunset are set to zero
    G0[!aman] <- D0[!aman] <- B0[!aman] <- kt[!aman] <- fd[!aman] <- 0

    result <- zoo(data.frame(kt, fd, G0, D0, B0), order.by=indSol)
    attr(result, 'match') <- mtch

    return(result)
}


