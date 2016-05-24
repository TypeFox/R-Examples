#####
decomp<-
function(Sigdata, ncomp=2, constant=TRUE, 
         typ=c("cw","lm"), control.args=list(),
         transf=TRUE, LEDpower=60, LEDwavelength=470,
         weight=FALSE, plot=TRUE, xylog=FALSE, lwd=3, outfile=NULL) {
    UseMethod("decomp")
} ###
### 2014.10.03.
decomp.default<-
function(Sigdata, ncomp=2, constant=TRUE, 
         typ=c("cw","lm"), control.args=list(),
         transf=TRUE, LEDpower=60, LEDwavelength=470,
         weight=FALSE, plot=TRUE, xylog=FALSE, lwd=3, outfile=NULL) {
    #### Stop if not.
    stopifnot(ncol(Sigdata)==2L,
              all(Sigdata[,2L,drop=TRUE]>0),
              length(ncomp)==1L, ncomp %in% seq(7L),
              length(constant)==1L, is.logical(constant),
              is.character(typ), all(typ %in% c("cw","lm")),
              class(control.args)=="list",
              all(names(control.args) %in% 
              list("factor","f","cr","maxiter","tol")),
              length(transf)==1L, transf==TRUE,
              length(LEDpower)==1L, is.numeric(LEDpower),
              length(LEDwavelength)==1L, is.numeric(LEDwavelength),
              length(weight)==1L, is.logical(weight),
              length(plot)==1L, is.logical(plot),
              length(xylog)==1L, is.logical(xylog),
              length(lwd)==1L, is.numeric(lwd),
              is.null(outfile) || is.character(outfile))
    ###
    tim<-as.numeric(Sigdata[,1L,drop=TRUE])
    sig<-as.numeric(Sigdata[,2L,drop=TRUE])
    ntim<-nrow(Sigdata)
    addc<-ifelse(constant==TRUE,1L,0L)
    n2<-2L*ncomp+addc
    uw<-ifelse(weight==FALSE,0L,1L)
    pars<-stdp<-vector(length=n2)
    typ1<-ifelse(typ[1L]=="cw",1L,2L)
    ### Check if data points is enough for fitting.
    if (ntim<n2) {
        stop("Error: data points is not enough for the model!")
    } # end if
    ### Default differential parameters.
    args<-list(factor=10L,f=0.5,cr=0.99,maxiter=500L,tol=0.1)
    args[names(control.args)]<-control.args
    factor<-args$factor
    f<-args$f
    cr<-args$cr
    maxiter<-args$maxiter
    tol<-args$tol
    stopifnot(is.numeric(factor), length(factor)==1L,
              factor>=5L, factor<=50L,
              is.numeric(f), length(f)==1L,
              f>0.0, f<=1.2,
              is.numeric(cr), length(cr)==1L,
              cr>0.0, cr<=1.0,
              is.numeric(maxiter), length(maxiter)==1L,
              maxiter>=10L, maxiter<=5000L,
              is.numeric(tol), length(tol)==1L, tol>0.0)
    ####
    fvec1<-vector(length=ntim)
    fmin<-0
    message<-0
    ###
    res<-.Fortran("decomp",as.double(tim),as.double(sig),as.integer(ntim), 
                  pars=as.double(pars),stdp=as.double(stdp),as.integer(n2),
                  as.integer(uw),as.integer(addc),as.integer(typ1),as.integer(factor), 
                  as.double(f),as.double(cr),as.integer(maxiter),as.double(tol), 
                  fvec1=as.double(fvec1),fmin=as.double(fmin),
                  message=as.integer(message),PACKAGE="numOSL")
    if (res$message!=0) {
        stop("Error: fail in decay curve decomposition!")
    } # end if
    ###
    h<-6.62606957e-34
    ny<-299792458/(LEDwavelength/10^9L)
    E<-h*ny
    LEDpower<-LEDpower/1000.0
    ###
    pars<-cbind(res$pars[seq(ncomp)], 
                res$stdp[seq(ncomp)],
                res$pars[(ncomp+1L):(2L*ncomp)], 
                res$stdp[(ncomp+1L):(2L*ncomp)],
                res$pars[(ncomp+1L):(2L*ncomp)]/LEDpower*E)
    pars<-pars[order(pars[,3L,drop=TRUE],decreasing=TRUE),,drop=FALSE]
    colnames(pars)<-c("Ithn", "Std.Ithn", "Lamda", "Std.Lamda", "Pcs")
    rownames(pars)<-paste("Comp.", seq(ncomp), sep="")
    ###
    if (constant==TRUE) {
        constant<-c("Constant"=res$pars[2L*ncomp+1L],
                    "Std.Constant"=res$stdp[2L*ncomp+1L])
    } else {
        constant<-0
    } # end if
    ###
    CompSig<-apply(cbind(pars[,1L,drop=TRUE], pars[,3L,drop=TRUE]), MARGIN=1L,
             function(x) if(typ[1L]=="cw") x[1L]*x[2L]*exp(-x[2L]*tim) else 
             x[1L]*x[2L]*(tim/max(tim))*exp(-x[2L]*tim^2L/2L/max(tim)))
    CompSig<-cbind(res$fvec1, CompSig)
    colnames(CompSig)<-c("Fit.Signal", paste("Comp.", seq(ncomp), sep=""))
    ### 
    if (plot==TRUE) {
        SigProp<-pars[,1L,drop=TRUE]/sum(pars[,1L,drop=TRUE])
        layout(cbind(c(1L,1L,1L,2L),c(1L,1L,1L,2L)))
        par(bg="grey95", mar=c(0,5.1,3.1,1.1))
        plot(tim, sig, log=ifelse(xylog==TRUE,"xy","x"), las=0, cex.main=1.5, lab=c(7,7,9),
             ylim=c(ifelse(xylog==TRUE,1e-5,0),ifelse(xylog==TRUE,max(sig)*2.0,max(sig)*1.1)), 
             cex.lab=1.5, xaxt="n", ylab="Photon Counts", xaxs="r", yaxs="i", type="p",
             pch=21, cex=ifelse(typ[1L]=="cw",1.5,1.25), bg="white", col="black")
        XaxisCentral<-median(axTicks(side=1L))
        colors<-c("blue", "green", "red", "deepskyblue", 
                  "purple", "orange", "brown")
        x<-seq(min(tim),max(tim),by=(max(tim)-min(tim))/ntim/100L)
        ###
        if (typ[1L]=="cw") {
             lines(x, eval(parse(text=paste("pars[",seq(ncomp),
                   ",1L,drop=TRUE]*pars[",seq(ncomp),
                   ",3L,drop=TRUE]*exp(-pars[",
                   seq(ncomp),",3L,drop=TRUE]*x)", 
                   collapse="+",sep="")))+constant[1L], 
                   lwd=lwd, col="black", lty="solid")
            if (constant[1L]>0) {
                points(tim, rep(constant[1L],ntim),
                       type="l", lty="dashed", lwd=lwd)
            } # end if
        } else if (typ[1L]=="lm") {
            lines(x, eval(parse(text=paste("pars[",seq(ncomp),
                  ",1L,drop=TRUE]*pars[", seq(ncomp),
                  ",3L,drop=TRUE]*(x/max(tim))*exp(-pars[",
                  seq(ncomp),",3L,drop=TRUE]*x^2/2L/max(tim))",
                  collapse="+",sep="")))+constant[1L]*x/max(tim), 
                  lwd=lwd, col="black", lty="solid")
            if (constant[1L]>0) {
                 points(tim, constant[1L]*tim/max(tim), 
                        type="l", lty="dashed", lwd=lwd)
            } # end if
        } # end if
        ###
        for (i in seq(ncomp)) {
            if (typ[1L]=="cw") {
                curve(pars[i,1L,drop=TRUE]*pars[i,3L,drop=TRUE]*
                      exp(-pars[i,3L,drop=TRUE]*x)+1e-7, 
                      lwd=lwd, col=colors[i], lty="solid", add=TRUE)
            } else if (typ[1L]=="lm") {
                curve(pars[i,1L,drop=TRUE]*pars[i,3L,drop=TRUE]*
                      (x/max(tim))*exp(-pars[i,3L,drop=TRUE]*x^2L/2L/max(tim))+1e-7, 
                      lwd=lwd, col=colors[i], lty="solid", add=TRUE)
            } # end if
        } # end for
        ###
        if (constant[1L]>0.0) {
            legend(ifelse(xylog==TRUE,"bottomleft",ifelse(typ[1L]=="cw",
                   "topright",ifelse(tim[which.max(sig)]>XaxisCentral, 
                   "topleft","topright"))),legend=c("Fitted.Curve",paste("Comp.", 
                   seq(ncomp)," (",round(SigProp*100,2L),"%)",sep=""),"Constant"),  
                   col=c("black",colors[seq(ncomp)],"black"),pch=c(21,rep(NA,ncomp+1L)), 
                   lty=c(rep("solid",ncomp+1L),"dashed"),yjust=2,ncol=1, 
                   cex=1.5,bty="o",lwd=lwd,pt.bg="white")
        } else {
            legend(ifelse(xylog==TRUE,"bottomleft",ifelse(typ[1L]=="cw",
                   "topright",ifelse(tim[which.max(sig)]>XaxisCentral, 
                   "topleft","topright"))),legend=c("Fitted.Curve",paste("Comp.", 
                   seq(ncomp)," (",round(SigProp*100,2L),"%)",sep="")),  
                   col=c("black",colors[seq(ncomp)]),pch=c(21,rep(NA,ncomp)), 
                   lty=c(rep("solid",ncomp+1L)),yjust=2,ncol=1, 
                   cex=1.5,bty="o",lwd=lwd,pt.bg="white")
        } # end if
        ### 
        grid(equilogs=FALSE)
        box(lwd=2L)
        ###
        par(mar=c(5.1,5.1,0,1.1))
        plot(tim, sig-CompSig[,1L,drop=TRUE], log="x", 
             las=0, lab=c(7,7,9), ylim=c(min(sig-CompSig[,1L,drop=TRUE])*1.05, 
             max(sig-CompSig[,1L,drop=TRUE])*1.05), xlab="Stimulated Time (s)", 
             cex.lab=1.5, ylab="Residuals", xaxs="r", yaxs="i", type="o", 
             pch=21, cex=ifelse(typ[1L]=="cw",0.75,0.5), bg="black", col="gray") 
        abline(h=0)
        box(lwd=2L)
        ###
        par(bg="transparent", mar=c(5,4,4,2)+0.1)
        layout(1L)
    } # end if
    ###
    output<-list("pars"=pars, 
                 "constant"=constant, 
                 "value"=res$fmin)
    ###
    if(!is.null(outfile)) {
        write.csv(CompSig, file=paste(outfile,".csv",sep=""))
    } # end if
    ###
    return(output)
} # end function decomp.default
#####
