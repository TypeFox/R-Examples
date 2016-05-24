checknPV <-
function(x)
{

namTAB<-rownames(x$outDAT)

SE<-x$outDAT[,"se"]
SP<-x$outDAT[,"sp"]
prev<-x$outDAT[,"prev"]
PPV0<-x$outDAT[,"PPV0"]
NPV0<-x$outDAT[,"NPV0"]

if(any(NPV0<=(1-prev))){
wNPV0w<-which(NPV0<=(1-prev))
ERRNPV0P<-paste("WARNING: The parameter(s) chosen in ",paste(namTAB[wNPV0w], collapse=", "), "contain(s) a threshold NPV0 which is smaller than 1-prevalence!")
warning(ERRNPV0P)
}

if(any(PPV0<=(prev))){
wPPV0w<-which(PPV0<=(prev))
ERRPPV0P<-paste("WARNING: The parameter(s) chosen in ", paste(namTAB[wPPV0w], collapse=", "), "contain(s) a threshold PPV0 which is smaller than the specified prevalence!") 
warning(ERRPPV0P)
}

NPVT<-x$outDAT[,"trueNPV"]
PPVT<-x$outDAT[,"truePPV"]

if(any(NPVT<=NPV0)){
whNPV<-which(NPVT<=NPV0)
ERRNPV<-paste("Hypothesized NPV is larger than the true NPV (resulting from your specifications of sensitivity, specificity and prevalence) in setting(s): ", paste(namTAB[whNPV], collapse=", "), ", i.e. H0 is true." , sep="")
warning(ERRNPV)
}

if(any(PPVT<=PPV0)){
whPPV<-which(PPVT<=PPV0)
ERRPPV<-paste("Hypothesized PPV is larger than the true PPV (resulting from your specifications of sensitivity, specificity and prevalence) in setting(s): ", paste(namTAB[whPPV], collapse=", "),
", i.e. H0 is true." , sep="")
warning(ERRPPV)
}

#MINN<-checkn.nPV(x=x)
#MAXLMIN<-nmax<MINN

#if(any(MAXLMIN))
#{
#wMAXLMIN<-which(MAXLMIN)
#WARNMAX<-paste("The specified maximal sample size (nmax) is too small to reject the specified hypothesis with specified power in setting (s) ",paste(namTAB[wMAXLMIN], collapse=", "), sep="")
#warning(WARNMAX)
#}

}

