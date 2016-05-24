change.pnparameters <-
structure(function 
                               (Asym = NA,
                               K = NA,
                                Infl = NA,
                                M = NA,
                                RAsym = NA,
                               Rk = NA,
                                Ri = NA,
                                RM = NA,
                                Amin = NA,
                                Amax = NA,
                                Kmin = NA,
                                Kmax = NA,
                                Imin = NA,
                                Imax = NA,
                                Mmin = NA,
                                Mmax = NA,
                                RAmin = NA,
                                RAmax = NA,
                                Rkmin = NA,
                                Rkmax = NA,
                                Rimin = NA,
                                Rimax = NA,
                                RMmin = NA,
                                RMmax = NA,                                
                                first.y = NA,
                                x.at.first.y = NA,
                                last.y = NA,
                                x.at.last.y = NA,
				twocomponent.x = NA,
				verbose = NA,
				force4par = NA,
				pn.options,
				Envir = .GlobalEnv
				) {
    newparams <- list(Asym = Asym, K = K, Infl = Infl, M = M, RAsym = RAsym,
        Rk = Rk, Ri = Ri, RM = RM, Amin = Amin, Amax = Amax, Kmin =Kmin, 
        Kmax = Kmax, Imin = Imin, Imax = Imax, Mmin = Mmin, Mmax = Mmax,
        RAmin = RAmin, RAmax = RAmax, Rkmin = Rkmin, Rkmax = Rkmax, 
        Rimin = Rimin, Rimax = Rimax, RMmin = RMmin, RMmax = RMmax, 
        first.y = first.y, x.at.first.y = x.at.first.y, 
        last.y = last.y,
        x.at.last.y = x.at.last.y, twocomponent.x = twocomponent.x, 
        verbose = verbose, force4par = force4par)
    pnoptnm<- as.character( pn.options[1] ) 
    if(as.character(list(Envir)) != "<environment>") stop ("No such environment")
    adjmodelparams <- try(get(pnoptnm, envir = Envir),silent=T)
    if( class(adjmodelparams) == "try-error" ) stop(paste(pnoptnm," not found in specified
    	environment. To search package environment do not specify an Envir value",sep=""))
     for (i in 1:9) {
        if (is.na(newparams[i]) == FALSE)
            adjmodelparams[names(adjmodelparams) %in% names(newparams[i])] <- newparams[i]
    }
     for (i in 25:31) {
            if (is.na(newparams[i]) == FALSE)
                adjmodelparams[names(adjmodelparams) %in% names(newparams[i])] <- newparams[i]
    }
    adjmodelparams <- adjmodelparams[names(adjmodelparams) %in% names(newparams)]
    pnmodelparams <- adjmodelparams 
    if(is.na(Amax)) Amax = pnmodelparams$Asym + (abs(pnmodelparams$Asym) * 2.5)
    if(is.na(Amin)) Amin = pnmodelparams$Asym - (abs(pnmodelparams$Asym) * 0.5)
    if(is.na(Kmax)) Kmax = pnmodelparams$K + (abs(pnmodelparams$K) * 0.5)
    if(is.na(Kmin)) Kmin = pnmodelparams$K - (abs(pnmodelparams$K) * 0.5)
    sImax <- Imax
    sImin <- Imin
    if(is.na(sImax)) Imax = pnmodelparams$Infl + (abs(pnmodelparams$Infl) * 2.5)
    if(is.na(sImin)) Imin = pnmodelparams$Infl - (abs(pnmodelparams$Infl) * 1.5)
    if(is.na(sImax))while (abs(Imax * Kmax) > 700) Imax = Imax * 0.9
    if(is.na(sImin))while (abs(Imin * Kmax) > 700) Imin = Imin * 0.9
    if(is.na(Mmax)) Mmax = pnmodelparams$M + abs(pnmodelparams$M * 2)
    if(is.na(Mmin)) Mmin = pnmodelparams$M - abs(pnmodelparams$M * 2)
    if (is.na(pnmodelparams$RAsym)) {
        pnmodelparams$RAsym <- pnmodelparams$Asym
        pnmodelparams$Rk <- pnmodelparams$K
        pnmodelparams$Ri <- pnmodelparams$Infl
        pnmodelparams$RM <- pnmodelparams$M
    }
    if(is.na(RAmax)) RAmax = pnmodelparams$RAsym + (abs(pnmodelparams$RAsym) *
        2.5)
    if(is.na(RAmin)) RAmin = pnmodelparams$RAsym - (abs(pnmodelparams$RAsym) *
        0.5)
    if(is.na(Rkmax)) Rkmax = pnmodelparams$Rk + (abs(pnmodelparams$Rk) * 0.5)
    if(is.na(Rkmin)) Rkmin = pnmodelparams$Rk - (abs(pnmodelparams$Rk) * 0.5)
    sRimax <- Rimax
    sRimin <- Rimin
    if(is.na(sRimax)) Rimax = pnmodelparams$Ri + (abs(pnmodelparams$Ri) * 1.25)
    if(is.na(sRimin)) Rimin = pnmodelparams$Ri - (abs(pnmodelparams$Ri) * 0.5)
    if(is.na(sRimax)) while (abs(Rimax * Rkmax) > 700) Rimax = Rimax * 0.9
    if(is.na(sRimax)) while (abs(Rimin * Rkmax) > 700) Rimin = Rimin * 0.9
    if(is.na(RMmax)) RMmax = pnmodelparams$RM + abs(pnmodelparams$RM * 2)
    if(is.na(RMmax)) RMmin = pnmodelparams$RM - abs(pnmodelparams$RM * 2)
    value2 <- c(Amin, Amax, Kmin, Kmax, Imin, Imax, Mmin, Mmax,
        RAmin, RAmax, Rkmin, Rkmax, Rimin, Rimax, RMmin, RMmax)
    skel1 <- rep(list(1), 16)
    value3 <- relist(value2, skel1)
    names(value3) <- c("Amin", "Amax", "Kmin", "Kmax", "Imin",
        "Imax", "Mmin", "Mmax", "RAmin", "RAmax", "Rkmin", "Rkmax",
        "Rimin", "Rimax", "RMmin", "RMmax")
    tmplate <- get(pnoptnm, envir = Envir)
    adjmodelparams[names(adjmodelparams) %in% names(value3)] <- value3[names(value3) %in% names(value3)]
    "%w/o%" <- function(x, y) x[!x %in% y] #--  x without y
    tmplate[names(tmplate) %in% names(adjmodelparams)] <- adjmodelparams[names(adjmodelparams) %in% names(adjmodelparams)]
    assign(pnoptnm, tmplate, Envir) 
    return(tmplate)
}
, ex = function(){
modpar(posneg.data$age,posneg.data$mass)
change.pnparameters(Asym=10000,Infl=80,M=5,RAsym=10000,Ri=240,RM=5)

change.pnparameters(M=1,RM=0.5,first.y=45.5)
}
)
