## we need ... because getMetric() will also pass parameters
## intended for other functions through ...
getTCP <-
function(x,
         TCPtcd50=NULL, TCPm=NULL, TCPn=NULL, TCPgamma50=NULL,
         EUDa=NULL, EUDfn=NULL, EUDab=NULL,
         TCPtype=c("probit", "logit", "poisson"), ...) {
    out <- getNTCP(x=x,
            NTCPtd50=TCPtcd50, NTCPm=TCPm, NTCPn=TCPn, NTCPgamma50=TCPgamma50,
            EUDa=EUDa, EUDfn=EUDfn, EUDab=EUDab,
            NTCPtype=TCPtype)
    outNames <- names(out)
    outNames[outNames == "NTCP"] <- "TCP"
    setNames(out, outNames)
}

getNTCP <-
function(x,
         NTCPtd50=NULL, NTCPm=NULL, NTCPn=NULL, NTCPgamma50=NULL,
         EUDa=NULL, EUDfn=NULL, EUDab=NULL,
         NTCPtype=c("probit", "logit", "poisson"), ...) {
    UseMethod("getNTCP")
}

getNTCP.DVHs <-
function(x,
         NTCPtd50=NULL, NTCPm=NULL, NTCPn=NULL, NTCPgamma50=NULL,
         EUDa=NULL, EUDfn=NULL, EUDab=NULL,
         NTCPtype=c("probit", "logit", "poisson"), ...) {
    NTCPtype <- match.arg(NTCPtype)
    stopifnot(!is.null(NTCPtd50),
              !is.null(NTCPm) || !is.null(NTCPgamma50),
              !is.null(NTCPn) || !is.null(EUDa))

    if(length(NTCPtd50) > 1) {
    	warning(paste0("Will only use NTCPtd50=", NTCPtd50[1]))
    	NTCPtd50 <- NTCPtd50[1]
    }

    if(NTCPtd50 <= 0) {
        warning("'NTCPtd50' must be > 0")
    	return(NA_real_)
    }

    if(is.null(NTCPm)) {
        if(length(NTCPgamma50) > 1) {
    		warning(paste0("Will only use NTCPgamma50=", NTCPgamma50[1]))
    		NTCPgamma50 <- NTCPgamma50[1]
        }

        if(NTCPgamma50 <= 0) {
            warning("'NTCPgamma50' must be > 0")
    		return(NA_real_)
        }

        NTCPm <- 1 / (NTCPgamma50*sqrt(2*pi))
    } else {
        if(length(NTCPm) > 1) {
    		warning(paste0("Will only use NTCPm=", NTCPm[1]))
    		NTCPm <- NTCPm[1]
        }

        if(NTCPm <= 0) {
            warning("'NTCPm' must be > 0")
    		return(NA_real_)
        }

        NTCPgamma50 <- 1 / (NTCPm*sqrt(2*pi))
    }

    if(!is.null(EUDa)) {
        if(length(EUDa) > 1) {
    		warning(paste0("Will only use EUDa=", EUDa[1]))
    		EUDa <- EUDa[1]
    	}

        if(isTRUE(all.equal(EUDa, 0))) {
    		warning("'EUDa' must not be zero")
    		return(NA_real_)
    	}

        NTCPn <- 1/EUDa
    } else {
        if(length(NTCPn) > 1) {
    		warning(paste0("Will only use NTCPn=", NTCPn[1]))
    		NTCPn <- NTCPn[1]
    	}

        if(is.infinite(NTCPn)) {
    		warning("'NTCPn' must not be infinite")
    		return(NA_real_)
    	}

        EUDa <- 1/NTCPn
    }

    EUD  <- getEUD(x, EUDa=EUDa, EUDfn=EUDfn, EUDab=EUDab)$EUD
    NTCP <- if(NTCPtype == "probit") {
        ## Lyman probit model based on EUD
        ## quantile at which to evaluate standard normal cdf
        qq <- (EUD-NTCPtd50) / (NTCPm*NTCPtd50)
        pnorm(qq, mean=0, sd=1, lower.tail=TRUE)
    } else if(NTCPtype == "logit") {
        ## Niemierko logit model based on EUD
        1 / (1 + ((NTCPtd50/EUD)^(4*NTCPgamma50)))
    } else if(NTCPtype == "poisson") {
        ## Kaellman Poisson model based on EUD
        2^(-exp(exp(1)*NTCPgamma50*(1-(EUD/NTCPtd50))))
    }

    data.frame(NTCP=NTCP,
               patID=x$patID,
               structure=x$structure,
               stringsAsFactors=FALSE)
}

getNTCP.DVHLst <-
function(x,
         NTCPtd50=NULL, NTCPm=NULL, NTCPn=NULL, NTCPgamma50=NULL,
         EUDa=NULL, EUDfn=NULL, EUDab=NULL,
         NTCPtype=c("probit", "logit", "poisson"), ...) {
    NTCPl <- Map(getNTCP, x,
                 NTCPtd50=list(NTCPtd50), NTCPm=list(NTCPm), NTCPn=list(NTCPn),
                 NTCPgamma50=list(NTCPgamma50), EUDa=list(EUDa),
                 EUDfn=list(EUDfn), EUDab=list(EUDab),
                 NTCPtype=list(NTCPtype))
    NTCPdf <- do.call("rbind", NTCPl)
    rownames(NTCPdf) <- NULL
    NTCPdf
}

getNTCP.DVHLstLst <-
function(x,
         NTCPtd50=NULL, NTCPm=NULL, NTCPn=NULL, NTCPgamma50=NULL,
         EUDa=NULL, EUDfn=NULL, EUDab=NULL,
         NTCPtype=c("probit", "logit", "poisson"), ...) {
    NTCPl <- Map(getNTCP, x,
                 NTCPtd50=list(NTCPtd50), NTCPm=list(NTCPm), NTCPn=list(NTCPn),
                 NTCPgamma50=list(NTCPgamma50), EUDa=list(EUDa),
                 EUDfn=list(EUDfn), EUDab=list(EUDab),
                 NTCPtype=list(NTCPtype))
    NTCPdf <- do.call("rbind", NTCPl)
    rownames(NTCPdf) <- NULL
    NTCPdf
}
