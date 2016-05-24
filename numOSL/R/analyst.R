#####
analyst<-
function(Sigdata, Redose, sig.channel=NULL, back.channel=NULL, 
         mr=0.01, typ="cw", nstart=100, upb=0.5, ErrorMethod=c("mc","sp"), 
         nsim=1000, weight=TRUE, plot=TRUE, model=NULL, origin=NULL) {
    UseMethod("analyst")
} #
### 2015.05.04, revised in 2016.01.20.
analyst.default<-
function(Sigdata, Redose, sig.channel=NULL, back.channel=NULL, 
         mr=0.01, typ="cw", nstart=100, upb=0.5, ErrorMethod=c("mc","sp"), 
         nsim=1000, weight=TRUE, plot=TRUE, model=NULL, origin=NULL) {
    ### Stop if not.
    stopifnot(nrow(Sigdata)>=25L, ncol(Sigdata)>=4L, ncol(Sigdata)%%2L==0L,
              is.vector(Redose), length(Redose)==(ncol(Sigdata)-2L)/2L,
              is.numeric(Redose), all(Redose>=0),
              is.null(sig.channel) || is.vector(sig.channel),
              all(sig.channel %in% seq(nrow(Sigdata))),
              is.null(back.channel) || is.vector(back.channel),
              all(back.channel %in% seq(nrow(Sigdata))),
              length(mr)==1L, is.numeric(mr), mr>=0, mr<=0.1,
              length(typ)==1L, is.character(typ), typ=="cw",
              is.numeric(nstart), length(nstart)==1L,
              nstart>=10L, nstart<=5000L,
              length(upb)==1L, is.numeric(upb), upb>0, upb<=10,
              is.character(ErrorMethod), all(ErrorMethod %in% c("mc","sp")),
              is.numeric(nsim), length(nsim)==1L, nsim>=100L, nsim<=3000L,
              length(weight)==1L, is.logical(weight),
              length(plot)==1L, is.logical(plot), 
              all(model %in% c("line","exp","lexp","dexp")),
              is.null(origin) || is.logical(origin),
              length(origin) %in% c(0L,1L))
    ###
    ndat<-nrow(Sigdata)
    Sigdata <- cbind(seq(ndat), Sigdata)
    ###
    if (is.null(sig.channel)) {
        sig.channel<-seq(4L)
    } # end if
    if (is.null(back.channel)) {
        back.channel<-(ndat-19L):ndat
    } # end if
    ###
    n<-length(sig.channel)
    m<-length(back.channel)
    k<-m/n
    ### 
    ### R function for calculating Lx/Tx.
    calLxTx<-function(sig1,sig2) {
        ###
        ### Total signal.
        Lx<-sum(sig1[sig.channel])
        Tx<-sum(sig2[sig.channel])
        ###
        ### Background signal.
        bLx<-n*mean(sig1[back.channel])
        bTx<-n*mean(sig2[back.channel])
        ###
        ### Net signal.
        netLx<-Lx-bLx      
        netTx<-Tx-bTx
        ### 
        if (abs(netLx)<=.Machine$double.eps)  netLx <- runif(n=1L, min=1e-8, max=1e-7)
        if (abs(netTx)<=.Machine$double.eps)  netTx <- runif(n=1L, min=1e-8, max=1e-7)
        ###
        ### Eqn.3 of Galbraith (2002).
        rsnetLx<-sqrt(Lx+bLx/k)/netLx
        rsnetTx<-sqrt(Tx+bTx/k)/netTx
        ###
        rsnetLx<-sqrt(rsnetLx^2L+mr^2L)
        rsnetTx<-sqrt(rsnetTx^2L+mr^2L)
        ###
        LxTx<-netLx/netTx
        sLxTx<-abs(LxTx)*sqrt(rsnetLx^2L+rsnetTx^2L)
        ###sLxTx<-abs(LxTx)*sqrt((Lx+bLx)/(Lx-bLx)^2L+(Tx+bTx)/(Tx-bTx)^2L)
        ###
        return(c(LxTx,sLxTx))
    } # end function calLxTx.
    ###
    ###
    nLxTx<-(ncol(Sigdata)-1L)/2L
    matLxTx<-matrix(nrow=nLxTx,ncol=2L)
    for (i in seq(nLxTx)) {
        matLxTx[i,]<-calLxTx(Sigdata[,2L*i,drop=TRUE],
                             Sigdata[,2L*i+1L,drop=TRUE])
    } # end for
    ###
    if (any(!is.finite(matLxTx))) {
        stop("Fail in calculating Lx/Tx!")
    } # end if
    ###
    Curvedata<-data.frame("Redose"=Redose,
                          "OSL"=matLxTx[-1L,1L,drop=TRUE],
                          "Std.OSL"=matLxTx[-1L,2L,drop=TRUE])
    NatureLxTx<-matLxTx[1L,,drop=TRUE]
    ###
    ### Calculate recycling ratio.
    lvl.dose<-as.numeric(levels(factor(Redose)))
    existrpd<-length(Redose)==length(lvl.dose)+1L
    if (existrpd==TRUE) {
        RepeatIndex<-apply(as.matrix(lvl.dose), MARGIN=1L, function(x,y) 
                     which(abs(x-y)<=.Machine$double.eps^0.5), Redose)
        RepeatIndex<-unlist(RepeatIndex[sapply(RepeatIndex,length)==2L])
        RecycleRatio<-Curvedata[,2L,drop=TRUE][RepeatIndex[2L]]/
                      Curvedata[,2L,drop=TRUE][RepeatIndex[1L]]
        #if (RecycleRatio>1.3 || RecycleRatio<0.7) {
            #cat("Note: recycling ratio is large (small)!\n")
        #} # end if
    } else {
        #cat("Note: recycling ratio is not available!\n")
        RecycleRatio<-NA
    } # end if
    ###
    ### Calculate recuperation.
    exist0d<-which(abs(Redose)<=.Machine$double.eps^0.5)
    if (length(exist0d)>0L) {
        Recuperation<-Curvedata[exist0d[1L],2L,drop=TRUE]/NatureLxTx[1L]
        #if (Recuperation>0.3) {
            #cat("Note: recuperation is large!\n")
        #} # end if
    } else {
        #cat("Note: recuperation is not available!\n")
        Recuperation<-NA
    } # end if  
    ###
    ###
    Models<-if(is.null(model)) c("line","exp") else model
    Origins<-if(is.null(origin)) c(TRUE,FALSE) else origin
    minf<-1e20
    ###
    for (i in Models) {
        for (j in Origins) {
            res<-try(calED(Curvedata=Curvedata, Ltx=NatureLxTx, model=i,
                           origin=j, nstart=50L, upb=upb, 
                           ErrorMethod=ErrorMethod, nsim=100L, 
                           weight=weight, plot=FALSE), silent=TRUE)
            ###
            if (class(res)!="try-error") {
                if (res$value<minf) {
                    model<-i
                    origin<-j
                    minf<-res$value
                    OK<-1L
                } # end if
            } # end if
        } # end for
    } # end for
    ###
    if (exists("OK")) {
        repeat {
            res<-try(calED(Curvedata=Curvedata, Ltx=NatureLxTx, model=model, 
                           origin=origin, nstart=nstart, upb=upb, 
                           ErrorMethod=ErrorMethod, nsim=nsim,
                           weight=weight, plot=plot), silent=TRUE)
            if (class(res)!="try-error") break
        } # end repeat.
    } else {
        stop("Error: fail in calculating ED!")
    } # end if   
    ###
    output<-list("Curvedata"=Curvedata, 
                 "Ltx"=NatureLxTx, 
                 "model"=paste(model,"(origin=",origin,")",sep=""), 
                 "LMpars"=res$LMpars, 
                 "value"=res$value, 
                 "ED"=res$ED, 
                 "RecyclingRatio"=RecycleRatio, 
                 "Recuperation"=Recuperation)
    ###
    return(output)
} # end function analyst.default.  
