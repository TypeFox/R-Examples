#####
fitGrowth<-
function(Curvedata, model=c("line","exp","lexp","dexp"),
         origin=FALSE, nstart=100, upb=0.5, 
         weight=TRUE, plot=TRUE) {
    UseMethod("fitGrowth")
} #
### 2014.10.02.
fitGrowth.default<-
function(Curvedata, model=c("line","exp","lexp","dexp"),
         origin=FALSE, nstart=100, upb=0.5, 
         weight=TRUE, plot=TRUE) {
    ### Stop if not.
    stopifnot(ncol(Curvedata)==3L,
              all(Curvedata[,1L,drop=TRUE]>=0),
              all(Curvedata[,3L,drop=TRUE]>0),
              is.character(model), length(model)>=1L,
              all(model %in% c("line","exp","lexp","dexp")),
              length(origin)==1L, is.logical(origin),
              is.numeric(nstart), length(nstart)==1L,
              nstart>=10L, nstart<=5000L,
              length(upb)==1L, is.numeric(upb), upb>0, upb<=10,
              length(weight)==1L, is.logical(weight),
              length(plot)==1L, is.logical(plot))
    ###
    dose<-as.numeric(Curvedata[,1L,drop=TRUE])
    doseltx<-as.numeric(Curvedata[,2L,drop=TRUE])
    sltx<-as.numeric(Curvedata[,3L,drop=TRUE])
    ndat<-nrow(Curvedata)
    n2<-if (model[1L]=="line") {
            1L+!origin
        } else if (model[1L]=="exp") {
            2L+!origin
        } else if (model[1L]=="lexp") {
            3L+!origin
        } else if (model[1L]=="dexp") {
            4L+!origin
    } # end if  
    ###
    if (ndat<n2) {
        stop("Error: data points is not enough for the model!")
    } # end if
    ###
    pars<-stdp<-vector(length=n2)
    model1<-if (model[1L]=="line") {
        0L } else if (model[1L]=="exp") {
        1L } else if (model[1L]=="lexp") {
        2L } else if (model[1L]=="dexp") {
        3L } # end if
    uw<-ifelse(weight==FALSE,0L,1L)
    fvec1<-vector(length=ndat)
    fmin<-0
    message<-0
    ###
    res<-.Fortran("fitGrowth",as.double(dose),as.double(doseltx),as.double(sltx),
                  as.integer(ndat),as.integer(n2),pars=as.double(pars),
                  stdp=as.double(stdp),as.double(upb),as.integer(model1),
                  as.integer(uw),as.integer(nstart),fvec1=as.double(fvec1),
                  fmin=as.double(fmin),message=as.integer(message),PACKAGE="numOSL")
    if (res$message!=0) {
        stop("Error: fail in growth curve fitting!")
    } # end if
    ###
    LMpars<-cbind(res$pars,res$stdp)
    colnames(LMpars)<-c("Pars","Std.Pars")
    rownames(LMpars)<-(c("a","b","c","d","e"))[seq(n2)]
    ###
    fit.value<-cbind(dose, doseltx, res$fvec1)
    colnames(fit.value)<-c("Redose", "Lx/Tx", "Fit.Lx/Tx")
    rownames(fit.value)<-paste("Redose", seq(ndat), sep="")
    ###
    output<-list("LMpars"=LMpars,
                 "value"=res$fmin,
                 "fit.value"=fit.value)
    ###
    if (plot==TRUE) {
        plot(dose, doseltx, main="Growth Curve", xlab="Dose (Gy)", 
             ylab="Standardised OSL", cex=2, cex.lab=1, cex.main=1.25) 
        x<-NULL
        if(origin==TRUE)  {
            if(model[1L]=="line") {
               curve(LMpars[1L,1L]*x, type="l", add=TRUE, lw=2)
            } else if(model[1L]=="exp") {
               curve(LMpars[1L,1L]*(1.0-exp(-LMpars[2L,1L]*x)),  
                     type="l", add=TRUE, lw=2)
            } else if(model[1L]=="lexp")  {
               curve(LMpars[1L,1L]*(1.0-exp(-LMpars[2L,1L]*x))+LMpars[3L,1L]*x, 
                     type="l", add=TRUE, lw=2)
            } else if(model[1L]=="dexp") {
               curve(LMpars[1L,1L]*(1.0-exp(-LMpars[2L,1L]*x))+
                     LMpars[3L,1L]*(1.0-exp(-LMpars[4L,1L]*x)), 
                     type="l", add=TRUE, lw=2)
            } # end if
        } else {
            if(model[1L]=="line") {
                curve(LMpars[1L,1L]*x+LMpars[2L,1L], type="l", 
                      add=TRUE, lw=2)
            } else if(model[1L]=="exp") {
                curve(LMpars[1L,1L]*(1.0-exp(-LMpars[2L,1L]*x))+LMpars[3L,1L], 
                      type="l", add=TRUE, lw=2)
            } else if(model[1L]=="lexp")  {
                curve(LMpars[1L,1L]*(1.0-exp(-LMpars[2L,1L]*x))+LMpars[3L,1L]*x+LMpars[4L,1L], 
                      type="l", add=TRUE, lw=2)
            } else if(model[1L]=="dexp") {
                curve(LMpars[1L,1L]*(1.0-exp(-LMpars[2L,1L]*x))+
                      LMpars[3L,1L]*(1.0-exp(-LMpars[4L,1L]*x))+LMpars[5L,1L], 
                      type="l", add=TRUE, lw=2)
            } # end if
        } # end if
        ###
        arrowsData<-Curvedata[Curvedata[,3L,drop=TRUE]>=1e-3,,drop=FALSE]
        options("warn"=-1)
        if (nrow(arrowsData)>=1L) {
            arrows(x0=arrowsData[,1L,drop=TRUE], 
                   y0=arrowsData[,2L,drop=TRUE]-arrowsData[,3L,drop=TRUE]/2L,
                   x1=arrowsData[,1L,drop=TRUE],
                   y1=arrowsData[,2L,drop=TRUE]+arrowsData[,3L,drop=TRUE]/2L,
                   code=3, lwd=2.5, angle=90, length=0.05, col="black")
        } # end if
        ###
        options("warn"=0)
        grid()
        box(lwd=2)
    } # end if
    return(output)
} # end function fitGrowth.default
#####
