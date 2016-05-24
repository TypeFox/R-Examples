# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
# ** Copyright UCAR (c) 1992 - 2004 
# ** University Corporation for Atmospheric Research(UCAR) 
# ** National Center for Atmospheric Research(NCAR) 
# ** Research Applications Program(RAP) 
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA 
# ** 2004/1/7 11:29:42 
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
table.stats<- function(obs, pred = NULL, fudge = 0.01, silent = FALSE) {
    ## internal function used in verify
    ## used with a binary forecast and a binary outcome.

    if(is.null(pred) & length(obs) ==4 ) {

        if(!silent){ print(" Assume data entered as c(n11, n01, n10, n00) Obs*Forecast") }
        a <- as.numeric( obs[1] )
        b <- as.numeric( obs[2] )
        c <- as.numeric( obs[3] )
        d <- as.numeric( obs[4] )
        tab.out <- matrix(c(a,c,b,d), nrow = 2)

    } 

    if(is.null(pred) & is.matrix(obs) & prod(dim(obs)) ==4 ) {

        if(!silent) print(" Assume contingency table has observed values in columns, forecasts in rows")
        obs <- as.numeric(obs)
        a <- obs[1] 
        b <- obs[3]
        c <- obs[2]
        d <- obs[4]
        tab.out <- matrix(c(a,c,b,d), nrow = 2)

    }
 
    if(!is.null(pred)& !is.null(obs)){
    
      tab.out <- table(as.numeric(obs), as.numeric(pred))
     
      a <- tryCatch(tab.out["1", "1"], error=function(e) 0) 
      b <- tryCatch(tab.out["0", "1"], error=function(e) 0) 
      c <- tryCatch(tab.out["1", "0"], error=function(e) 0)
      d <- tryCatch(tab.out["0", "0"], error=function(e) 0)
      # a <-tab.out["1","1"]
      # b <-tab.out["0","1"]
      # c <-tab.out["1","0"]
      # d <-tab.out["0","0"]

    } ## close else
    
    
    ###
    n <- a + b + c + d
    if(n == 0) n <- fudge

    s <- (a+c)/n
    
    TS  <- a /(a+b+c + fudge)
    POD<- H <- a /(a+c + fudge)  ## hit rate
    F   <- b /(b+d + fudge)  ## false alarm rate
    TS.se <- sqrt((TS^2)*((1-H)/(a + fudge) + b*(1-F)/((a+b+c)^2 + fudge)))
    
    SH2 <- H*(1-H)/(a+c + fudge)
    SF2 <- F*(1-F)/(b+d + fudge)
    POD.se <- sqrt(SH2)
    F.se <- sqrt(SF2)
    
    M   <- c /(a+c + fudge)  ## miss rate
    FAR <- b/(a+b + fudge)  ## false alarm ratio
    FAR.se <- sqrt((FAR^4)*((1-H)/(a + fudge) + (1-F)/(b + fudge))*(a^2)/(b^2 + fudge))
    HSS <- 2*(a*d - b*c)/ (1.0*(a+c)*(c+d) + 1.0*(a+b)*(b+d) + fudge)  
    SHSS2 <- SF2*(HSS^2)*(1/(H-F + fudge) + (1-s)*(1-2*s))^2 + SH2*(HSS^2)*(1/(H-F + fudge) - s*(1-2*s))^2
    HSS.se = sqrt(SHSS2)
    PSS <- 1 - M  - F ## Pierce skill score
    PSS.se <- sqrt(SH2 + SF2)
    KSS <- (a*d - b*c)/((a+c)*(b + d) + fudge) ## similar to Pierc
    PC <- (a+d)/(a+b+c+d + fudge)
    PC.se <- sqrt(s*H*(1-H)/n + (1-s)*F*(1-F)/n)

    if(a + c == 0) BIAS <- (a+b)/fudge
    else BIAS <- (a + b)/(a + c)

    if(b * c == 0) OR   <- a*d/fudge ## odds ratio
    else OR <- a * d / (b * c)

    if(a * b + b * c == 0) ORSS <- (a*d - b*c)/ fudge
    else ORSS <- (a*d - b*c)/ (a*d + b*c ) ## odds ratio skill score

    HITSrandom <- 1.0* (a+c)*(a+b)/n
    p <- (a+c)/n
   
    if(a + b + c - HITSrandom == 0) ETS <- (a - HITSrandom)/fudge 
    else ETS <- (a-HITSrandom)/(a+b+c-HITSrandom)

    if(2 - HSS == 0) ETS.se <- sqrt(4 * SHSS2/fudge)
    else ETS.se <- sqrt(4*SHSS2/((2-HSS)^4))

    if(b * c == 0) theta <- a * d / fudge
    else theta <- (a*d)/(b*c)

    log.theta <- log(a) + log(d) - log(b) - log(c) 

    if(a == 0) a.z <- fudge
    else a.z <- a

    if(b == 0) b.z <- fudge
    else b.z <- b

    if(c == 0) c.z <- fudge
    else c.z <- c

    if(d == 0) d.z <- fudge
    else d.z <- d

    if(1/a.z + 1/b.z + 1/c.z + 1/d.z == 0) n.h <- 1/fudge
    else n.h <- 1/( 1/a.z + 1/b.z + 1/c.z + 1/d.z)

    if(theta + 1 == 0) yules.q <- (theta - 1)/fudge
    else yules.q <- (theta - 1)/(theta + 1) 

    if(n.h == 0) SLOR2 <- 1/fudge
    else SLOR2 <- 1/n.h

    LOR.se <- sqrt(SLOR2)

    if(OR + 1 == 0) ORSS.se <- sqrt(SLOR2 * 4 * OR^2 / fudge)
    else ORSS.se <- sqrt(SLOR2*4*OR^2/((OR+1)^4))

    if(log(a/n) == 0) {

	eds <- 2 * log((a + c)/n)/fudge - 1
	seds <- (log((a+b)/n)+log((a+c)/n)) / fudge - 1

    } else {

	eds <- 2*log((a+c)/n)/log(a/n) - 1
	seds <- (log((a+b)/n)+log((a+c)/n)) /log(a/n) - 1

    }

    eds.se <- 2*abs(log(p))/(H*(log(p) + log(H))^2)*sqrt(H*(1-H)/(p*n))
    seds.se <- sqrt(H*(1-H)/(n*p)) *(-log(BIAS*p^2)/(H*log(H*p)^2))

    if(log(F) + log(H) == 0) EDI <- (log(F) - log(H)) / fudge
    else EDI <- (log(F) - log(H))/(log(F) + log(H))

    EDI.se <- 2*abs(log(F) + H/(1-H)*log(H))/(H*(log(F) + log(H))^2)*sqrt(H*(1-H)/(p*n))
    SEDI <- (log(F) - log(H) - log(1-F) + log(1-H))/(log(F) + log(H) + log(1-F) + log(1-H))
    SEDI.se <- 2*abs(((1-H)*(1-F)+H*F)/((1-H)*(1-F))*log(F*(1-H)) + 2*H/(1-H)*log(H*(1-F)))/(H*(log(F*(1-H)) + log(H*(1-F)))^2)*sqrt(H*(1-H)/(p*n))
    
    return(list(tab = tab.out, TS = TS, TS.se = TS.se,
                POD = POD, POD.se = POD.se, M = M,  F = F,  F.se = F.se, FAR = FAR , FAR.se = FAR.se, HSS = HSS, HSS.se = HSS.se,
                PSS = PSS, PSS.se = PSS.se, KSS = KSS,
                PC = PC, PC.se = PC.se,  BIAS = BIAS, ETS = ETS, ETS.se = ETS.se, theta = theta, log.theta = log.theta, LOR.se = LOR.se, n.h = n.h, orss = yules.q, orss.se = ORSS.se,
    	    eds = eds, eds.se=eds.se, seds = seds, seds.se = seds.se,
    	    EDI = EDI, EDI.se = EDI.se, SEDI = SEDI, SEDI.se = SEDI.se))
    
} # end of function.


