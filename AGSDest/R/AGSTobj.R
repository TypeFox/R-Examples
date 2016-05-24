#' The \code{AGSTobj} includes design and outcome of primary and secondary trial.
#'
#' A \code{AGSTobj} object is designed.
#'
#' The function \code{summary} returns an object of \code{class} \code{AGSTobj}.
#' \code{ctype} defines the type of confidence interval that is calculated.
#' \tabular{ll}{
#' \code{"r"} \tab Repeated confidence bound for a GSD with design adaptations\cr
#' \code{"so"} \tab Confidence bound for a GSD with design adaptation based on the stage-wise ordering\cr
#' }
#' The calculated confidence bounds are saved as:
#' \tabular{ll}{
#' \code{cb.r} \tab repeated confidence bound \cr
#' \code{cb.so} \tab confidence bound based on the stage-wise ordering \cr
#' }
#' \code{ptype} defines the type of p-value that is calculated.
#' \tabular{ll}{
#' \code{"r"}\tab Repeated p-value for a GSD with design adaptations \cr
#' \code{"so"}\tab Stage-wise adjusted p-value for a GSD with design adaptations \cr' 
#' }
#' The calculated p-values are saved as:
#' \tabular{ll}{
#' \code{pvalue.r}\tab repeated p-value \cr
#' \code{pvalue.so}\tab stage-wise adjusted p-value\cr
#' }
#' \code{etype} defines the type of point estimate
#'    \tabular{ll}{
#'    \code{"ml"}\tab maximum likelihood estimate (ignoring the sequential and adaptive nature of the design)\cr
#'    \code{"mu"}\tab median unbiased estimate (stage-wise lower confidence bound at level 0.5) for a GSD with design adaptations\cr
#'    \code{"cons"}\tab conservative estimate (repeated lower confidence bound at level 0.5) for a GSD with design adaptations\cr
#'    }
#'  
#'  The calculated point estimates are saved as:
#'  \tabular{ll}{
#'  \code{est.ml}\tab Maximum likelihood estimate\cr
#'  \code{est.mu}\tab Median unbiased estimate\cr
#'  \code{est.cons}\tab Conservative estimate\cr
#'  }
#'  
#'  The stage-wise adjusted confidence bound, p-value and the median unbiased point estimate can only be calculated at the stage where the trial stops and are only valid if the stopping rule is met.
#'
#'  The repeated confidence bound, repeated p-value, conservative estimate and maximum likelihood estimate can be calculated at every stage of the trial and
#'  not just at the stage where the trial stops and are also valid if the stopping rule is not met. 
#'  For calculating the repeated confidence bounds or p-values the user has to specify \code{sTo} (secondary trial outcome) in the object \code{AGSTobj} (see example below). If the stopping rule is not met in object \code{sTo} then sta e-wise adjusted confidence bounds and p-values will not be computed while a warning message is given when their computation have erroneously been specified.
#' 
#' @title Adaptive group sequential trial object (AGSTobj)
#' @param x object of the \code{class} \code{AGSTobj}
#' @param object object of the \code{class} \code{AGSTobj}
#' @param main Title of the plots (default: first plot: "primary trial"; second plot: "secondary trial")
#' @param print.pdf option; if TRUE a pdf file is created. Instead of setting print.pdf to TRUE, the user can specify a character string giving the name or the path of the file.
#' @param ctype confidence type: repeated "r" or stage-wise ordering "so" (default: c("r", "so"))
#' @param ptype p-value type: repeated "r" or stage-wise ordering "so" (default: c("r", "so"))
#' @param etype point estimate: maximum likelihood "ml", median unbiased "mu", or conservative "cons" (default: c("ml", "mu", "cons"))
#' @param overwrite option; if TRUE all old values are deleted and new values are calculated (default: FALSE)
#' @param ... additional arguments.
#' @return An object of class \code{AGSTobj}, which is basically a list with the elements
#' \item{cb.so}{ confidence bound based on the stage-wise ordering (stage-wise adjusted confidence bound)}
#' \item{cb.r}{ repeated confidence bound}
#' \item{pvalue.so}{ p-value based on the stage-wise ordering (stage-wise adjusted p-value)}
#' \item{pvalue.r}{ repeated p-value}
#' \item{est.ml}{ maximum likelihood estimate}
#' \item{est.mu}{ median unbiased point estimate}
#' \item{est.cons}{ conservative point estimate}
#' \item{pT}{}
#' \item{K}{ number of stages}
#' \item{al}{ alpha (type I error rate)}
#' \item{a}{ lower critical bounds of primary group sequential design (are currently always set to -8)}
#' \item{b}{ upper critical bounds of primary group sequential design}
#' \item{t}{ vector with cumulative information fraction}
#' \item{SF}{ spending function (for details see below)}
#' \item{phi}{ parameter of spending function when SF=3 or 4 (for details see below)}
#' \item{alab}{ alpha-absorbing parameter values of primary group sequential design}
#' \item{als}{ alpha-values ''spent'' at each stage of primary group sequential design}
#' \item{Imax}{ maximum information number}
#' \item{delta}{ effect size used for planning the primary trial}
#' \item{cp}{ conditional power for planning the primary trial}
#' \item{iD}{}
#' \item{L}{ stage of the adaptation}
#' \item{z}{ z-statistic at adaptive interim analysis}
#' \item{sT}{}
#' \item{K}{ number of stages }
#' \item{al}{ conditional rejection probability}
#' \item{a}{ lower critical bounds of secondary group sequential design (are currently always set to -8)}
#' \item{b}{ upper critical bounds of secondary group sequential design}
#' \item{t}{ vector with cumulative information fraction}
#' \item{SF}{ spending function (for details see below)}
#' \item{phi}{ parameter of spending function when SF=3 or 4 (for details see below)}
#' \item{Imax}{ maximum information number}
#' \item{delta}{ effect size used for planning the secondary trial}
#' \item{cp}{ conditional power for planning the secondary trial}
#' \item{sTo}{}
#' \item{T}{ stage where trial stops }
#' \item{z}{ z-statistic at stage where trial stops}
#'
#' @author Niklas Hack \email{niklas.hack@@meduniwien.ac.at} and Werner Brannath \email{werner.brannath@@meduniwien.ac.at}
#' @note
#' The \code{AGSTobj} should always have the same ordering and names as given in the list above or as given in the example.
#' 1. \code{pt}, 2. \code{iD}, 3. \code{sT}, 4. \code{sTo}
#' \code{SF} defines the spending function.
#' \code{SF} = 1 O'Brien and Fleming type spending function of Lan and DeMets (1983)
#' \code{SF} = 2 Pocock type spending function of Lan and DeMets (1983)
#' \code{SF} = 3 Power family (\eqn{c_\alpha* t^\phi}); phi must be greater than 0.
#' \code{SF} = 4 Hwang-Shih-DeCani family;\eqn{(1-e^{-\phi t})/(1-e^{-\phi})}, where \code{phi} cannot be 0.
#' A value of \code{SF}=3 corresponds to the power family. Here, the spending function is \eqn{t^{\phi}},
#' where phi must be greater than 0. A value of \code{SF}=4 corresponds to the Hwang-Shih-DeCani family,
#' with the spending function \eqn{(1-e^{-\phi t})/(1-e^{-\phi})}, where phi cannot be 0.
#' If a path is specified for \code{print.pdf}, all \\ must be changed to /. If a filename is specified the ending of the file must be (.pdf).
#' In the current version the vector of lower bounds \code{a} should be set to rep(-8,K)
#' 
#' @seealso \code{\link{AGSTobj}}, \code{\link{print.AGSTobj}}, \code{\link{plot.AGSTobj}}, \code{\link{summary.AGSTobj}}
#' @examples
#' \dontrun{
#' pT=plan.GST(K=3,SF=4,phi=-4,alpha=0.05,delta=6,pow=0.9,compute.alab=TRUE,compute.als=TRUE)
#' 
#' iD=list(T=1, z=1.090728)
#' 
#' swImax=0.0625
#' 
#' I2min=3*swImax
#' I2max=3*swImax
#'
#' sT=adapt(pT=pT,iD=iD,SF=1,phi=0,cp=0.8,theta=5,I2min,I2max,swImax)
#'
#' sTo=list(T=2, z=2.393)
#'
#' AGST<-as.AGST(pT=pT,iD=iD,sT=sT,sTo=sTo)
#' AGST
#' plot(AGST)
#'
#' AGST<-summary(AGST)
#' plot(AGST)
#'
#' ##The repeated confidence interval and p-value at an earlier stage
#' ##than the one where the trial stops (T=3).
#'
#' summary(as.AGST(pT,iD,sT,sTo=list(T=1,z=1.7)),ctype="r",ptype="r")
#'
#' ##If the stage-wise adjusted confidence interval is calculated at this stage,
#' ##the function returns an error message
#' summary(as.AGST(pT,iD,sT,sTo=list(T=1,z=1.7)),ctype="so",ptype="so")
#' }
#' 
#' @keywords datasets 
AGSTobj <- function(x, ...) UseMethod("AGSTobj")

#' @rdname AGSTobj
#' @export
plot.AGSTobj<-function(x, main=c("primary trial", "secondary trial"), print.pdf=FALSE, ...) {

    if(!is.null(x[[1]]$t) && !is.null(x[[2]]$T) && !is.null(x[[3]]$t)) {
                                        #par.old=par(no.readonly = TRUE)
                                        #par(mfrow=c(2,1),oma = c(0,12,0,0))




        if(mode(print.pdf) == "logical" && print.pdf == FALSE) {
            nf <- layout(matrix(c(1,2,3,4),2,2, byrow=TRUE), c(2,3), c(3,3))
            layout.show(nf)
        }
        if(mode(print.pdf) == "logical" && print.pdf == TRUE) pdf("AGST.pdf")
        if(mode(print.pdf) == "character") pdf(print.pdf)

        plot(1,1, type="n", ann=FALSE, axes=FALSE, col="black")
        box(col="black")

        text(1, 1, paste(length(x[[1]]$t)," stage primary trial,\n",switch(x[[1]]$SF,"OBF","Pocock","Power family","Hwang-Shih-DeCani","self specified","NA")," boundaries\n at level ",x[[1]]$al,
                         if(!is.null(x[[1]]$Imax))paste(",\nImax=",round(x[[1]]$Imax,digits=2)),
                         if(!is.null(x[[1]]$delta))if(x[[1]]$delta!=0)paste(",\ndelta=",round(x[[1]]$delta,digits=2)),
                         if(!is.null(x[[2]]$T))paste("\n\nAdaptation at stage T=",x[[2]]$T), " \nwith z = ",round(x[[2]]$z,digits=3)),cex=0.8)
        
        plot(x[[1]]$t,x[[1]]$b,axes=FALSE,ylim=c(0,max(x[[1]]$b)+1.5),xlab="Cumulative Information Fraction",ylab="Wald Teststatistic",cex.lab=1.25,type="l",
             main=main[1], lwd=1.5,cex.main=0.9,cex.lab=0.9,cex=0.9,cex.axis=0.9,cex.sub=0.9,...)

        axis(2,cex.axis=0.9)
        axis(1, at=x[[1]]$t, labels=round(x[[1]]$t,digits=3),cex.axis=0.9)
        points(x[[1]]$t,x[[1]]$b,lt=1,pch=19)
        points(x[[1]]$t[x[[2]]$T],x[[2]]$z,pch=21)

                                        #if(print.pdf){dev.off();cat("\nFile pT.pdf is created");pdf("sT.pdf")}
                                        #else{print("Click to switch to next plot");par(ask=TRUE)}

                                        #if(locator(1)){

        plot(1,1, type="n",ann=FALSE, axes=FALSE, col="black")
        box(col="black")

        text(1,1,paste(length(x[[3]]$t)," stage secondary trial,\n",switch(x[[3]]$SF,"OBF","Pocock","Power family","Hwang-Shih-DeCani","self specified","NA")," boundaries\n at level ",round(x[[3]]$al,digits=3),
                       if(!is.null(x[[3]]$Imax))paste(",\nImax=",round(x[[3]]$Imax,digits=2)),
                       if(!is.null(x[[3]]$delta))if(x[[3]]$delta!=0)paste(",\ndelta=",round(x[[3]]$delta,digits=2)),
                       if(!is.null(x[[4]]$T) && x[[4]]$T!=0){paste("\n\nTrial stops at stage T=",x[[4]]$T,"\nwith z = ",round(x[[4]]$z,digits=3))},
                       paste("\n"),
                                        #if(!is.null(x$pvalue.r) || !is.null(x$pvalue.so))paste("\n"),
                       if(!is.null(x$pvalue.r))paste("\npvalue.r = ",round(x$pvalue.r,digits=3)),
                       if(!is.null(x$pvalue.so))paste("\npvalue.so = ",round(x$pvalue.so,digits=3)),
                                        #if(!is.null(x$cb.r) || !is.null(x$cb.so))paste("\n"),
                       if(!is.null(x$cb.r))paste("\ncb.r = ",round(x$cb.r,digits=3)),
                       if(!is.null(x$cb.so))paste("\ncb.so = ",round(x$cb.so,digits=3)),  
                                        #if(!is.null(x$est.ml) || !is.null(x$est.mu))paste("\n"),
                       if(!is.null(x$est.ml))paste("\nest.ml = ",round(x$est.ml,digits=3)),
                       if(!is.null(x$est.mu))paste("\nest.mu = ",round(x$est.mu,digits=3))),cex=0.8)


        plot(x[[3]]$t,x[[3]]$b,axes=FALSE,ylim=c(0,max(x[[3]]$b)+1.5),xlab="Cumulative Information Fraction",ylab="Wald Teststatistic",cex.lab=1.25,type="l",
             main=main[2], lwd=1.5,cex.main=0.9,cex.lab=0.9,cex=0.9,cex.axis=0.9,cex.sub=0.9,...)

        axis(2,cex.axis=0.9)
        axis(1, at=x[[3]]$t, labels=round(x[[3]]$t,digits=3),cex.axis=0.9)
        points(x[[3]]$t,x[[3]]$b,lt=1,pch=19)

        if(!is.null(x[[4]]$T) && x[[4]]$T!=0)points(x[[3]]$t[x[[4]]$T],x[[4]]$z,pch=21)

        if(mode(print.pdf) == "logical" && print.pdf == TRUE) {
            dev.off()
            cat("File AGST.pdf is created\n")
        }
        if(mode(print.pdf) == "character") {
            dev.off()
            cat("File", print.pdf, "is created\n")
        }
    }
}

#' @rdname AGSTobj
#' @export
print.AGSTobj <- function(x, ...) {
    
    if(!is.null(x$cb.r)) cat(paste("repeated lower confidence bound: ",round(x$cb.r,digits=3),"\n\n"))
    if(!is.null(x$cb.so)) cat(paste("stage-wise adjusted lower confidence bound: ",round(x$cb.so,digits=3),"\n\n"))
   
    if(!is.null(x$pvalue.r)) cat(paste("repeated p-value: ",round(x$pvalue.r,digits=3),"\n\n"))
    if(!is.null(x$pvalue.so)) cat(paste("stage-wise adjusted p-value: ",round(x$pvalue.so,digits=3),"\n\n"))
    
    if(!is.null(x$est.ml)) cat(paste("maximum likelihood estimate: ",round(x$est.ml,digits=3),"\n\n"))
    if(!is.null(x$est.mu)) cat(paste("median unbiased estimate: ",round(x$est.mu,digits=3),"\n\n"))
    if(!is.null(x$est.cons)) cat(paste("conservative estimate: ",round(x$est.cons,digits=3),"\n\n"))
    
    cat(paste("Primary trial: \n\n"))
    cat(paste(length(x[[1]]$t)," stage group sequential design"))
    cat(paste("\n",expression(alpha),": ",round(x[[1]]$al,digits=3),"  SF: ",x[[1]]$SF,"  phi: ",x[[1]]$phi,if(!is.null(x[[1]]$Imax))paste("  Imax: ",round(x[[1]]$Imax,digits=2)),if(!is.null(x[[1]]$delta))if(x[[1]]$delta!=0)paste("  delta: ",round(x[[1]]$delta,digits=2)),if(!is.null(x[[1]]$cp))if(x[[1]]$cp!=0)paste("  cp: ",round(x[[1]]$cp,digits=2))))

    cat("\n")
    mat <- matrix(c(round(x[[1]]$b,digits=3),round(x[[1]]$a,digits=3),round(x[[1]]$t,digits=3)),ncol=length(x[[1]]$t),byrow=TRUE)
    colnames(mat) <- rep("",length(x[[1]]$t))
    rownames(mat) <- c("Upper bounds","Lower bounds","Information fraction")
    print(mat)
    
    cat("\n");
    if(!is.null(x[[1]]$alab)) {
        mat_als <- matrix(c(round(x[[1]]$als,digits=3),round(x[[1]]$alab[1:length(x[[1]]$t)],digits=3)),ncol=length(x[[1]]$t),byrow=TRUE)
        colnames(mat_als) <- rep("",length(x[[1]]$t))
        rownames(mat_als) <- c("als","alab")
        print(mat_als)
        cat("\n")
    }
    
    cat("\n\n")
    cat(paste("interim data: \n"))
    cat(paste("\n\tT: ", x[[2]]$T, "  z: ", round(x[[2]]$z, digits=3)))
    
    cat(paste("\n\nSecondary trial: \n\n"))
    cat(paste(length(x[[3]]$t)," stage group sequential design"))
    cat(paste("\n cer: ",round(x[[3]]$al,digits=3),"  SF: ",x[[3]]$SF,"  phi: ",x[[3]]$phi,if(!is.null(x[[3]]$Imax))paste("  Imax: ",round(x[[3]]$Imax,digits=2)),if(!is.null(x[[3]]$delta))if(x[[3]]$delta!=0)paste("  delta: ",round(x[[3]]$delta,digits=2)),if(!is.null(x[[3]]$cp))if(x[[3]]$cp!=0)paste("  cp: ",round(x[[3]]$cp,digits=2))))

    
    cat("\n")
    mat <- matrix(c(round(x[[3]]$b,digits=3),round(x[[3]]$a,digits=3),round(x[[3]]$t,digits=3)),ncol=length(x[[3]]$t),byrow=TRUE)
    colnames(mat) <- rep("",length(x[[3]]$t))
    rownames(mat) <- c("Upper bounds","Lower bounds","Information fraction")
    print(mat)
    
    if(!is.null(x[[4]]$T) && x[[4]]$T != 0) {
        cat(paste("\n\nSecondary trial outcome: \n"))
        cat(paste("\n\tT: ", x[[4]]$T, "  z: ", round(x[[4]]$z, digits=3), "\n\n"))
    }
    else cat(paste("\n"))
}

#' @rdname AGSTobj
#' @export
summary.AGSTobj <- function(object, ctype=c("r", "so"), ptype=c("r", "so"), etype=c("ml", "mu", "cons"), overwrite=FALSE, ...) {

    ## repeated confidence bound
    if("r" %in% ctype) {
        if(is.null(object$cb.r) || object$cb.r==0 || (!is.null(object$cb.r) && overwrite)) {
            object$cb.r <- cb.r.ad(object$pT, object$iD, object$sT, object$sTo)
        }    
    }

    ## stage-wise adjusted confidence bound
    if("so" %in% ctype) {
        if(is.null(object$cb.so) || object$cb.so==0 || (!is.null(object$cb.so) && overwrite)) {
            if(object$sTo$z < object$sT$b[object$sTo$T] && object$sTo$T < object$sT$K) cat("cb.so : z < b[T]; Stopping rule NOT met.\n")
            else object$cb.so <- cb.so.ad(object$pT, object$iD, object$sT, object$sTo)
        }
    }

    ## repeated p-value
    if("r" %in% ptype) {
        if(is.null(object$pvalue.r) || object$pvalue.r==0 || (!is.null(object$pvalue.r) && overwrite)) {
            object$pvalue.r <- P.r.ad(h=0, object$pT, object$iD, object$sT, object$sTo)
        }
    }

    ## stage-wise adjusted p-value
    if("so" %in% ptype) {
        if(is.null(object$pvalue.so) || object$pvalue.so==0 || (!is.null(object$pvalue.so) && overwrite)) {
            if(object$sTo$z < object$sT$b[object$sTo$T]) cat("pvalue.so : z < b[T]; Stopping rule NOT met.\n")
            else object$pvalue.so <- P.so.ad(h=0, object$pT, object$iD, object$sT, object$sTo)
        }

        ## maximum likelihood estimate
        if("ml" %in% etype) {
            if(is.null(object$est.ml) || object$est.ml==0  || (!is.null(object$est.ml) && overwrite)) {
                object$est.ml <- (object$iD$z*sqrt(object$pT$t[object$iD$T]*object$pT$Imax) + object$sTo$z*sqrt(object$sT$t[object$sTo$T]*object$sT$Imax)) /
                    (object$pT$t[object$iD$T]*object$pT$Imax+object$sT$t[object$sTo$T]*object$sT$Imax)
            }
        }

        ## median unbiased estimate
        if("mu" %in% etype) {
            if(is.null(object$est.mu) || object$est.mu==0 || (!is.null(object$est.mu) && overwrite)) {
                if(object$sTo$z < object$sT$b[object$sTo$T] && object$sTo$T < object$sT$K) cat("est.mu : z < b[T]; Stopping rule NOT met.\n")
                else object$est.mu <- cb.so.ad(object$pT, object$iD, object$sT, object$sTo, level=0.5)
            }
        }

        ## conservative estimate
        if("cons" %in% etype) {
            if(is.null(object$est.cons) || object$est.cons==0 || (!is.null(object$est.cons) && overwrite)) {
                object$est.cons <- cb.r.ad(object$pT, object$iD, object$sT, object$sTo, level=0.5)
            }
        }
    }
    
    object
}

#' @rdname AGSTobj
#' @export
print.summary.AGSTobj <- function(x,...) {
        
    if(!is.null(x$type)) cat(paste("type: ",x$type,"\n\n"))        
    if(!is.null(x$cb.r)) cat(paste("repeated lower confidence bound: ",round(x$cb.r,digits=3),"\n\n"))
    if(!is.null(x$cb.so)) cat(paste("stage-wise adjusted lower confidence bound: ",round(x$cb.so,digits=3),"\n\n"))
    
    if(!is.null(x$pvalue.r)) cat(paste("repeated p-value: ",round(x$pvalue.r,digits=3),"\n\n"))
    if(!is.null(x$pvalue.so)) cat(paste("stage-wise adjusted p-value: ",round(x$pvalue.so,digits=3),"\n\n"))
    
    if(!is.null(x$est.ml)) cat(paste("maximum likelihood estimate: ",round(x$est.ml,digits=3),"\n\n"))
    if(!is.null(x$est.mu)) cat(paste("median unbiased estimate: ",round(x$est.mu,digits=3),"\n\n"))
    if(!is.null(x$est.cons)) cat(paste("conservative estimate: ",round(x$est.cons,digits=3),"\n\n"))
    
    cat(paste("Primary trial: \n\n"))
    cat(paste(length(x[[1]]$t)," stage group sequential design"))
    cat(paste("\n",expression(alpha),": ",round(x[[1]]$al,digits=3),"  SF: ",x[[1]]$SF,"  phi: ",x[[1]]$phi,if(!is.null(x[[1]]$Imax))paste("  Imax: ",round(x[[1]]$Imax,digits=2)),if(!is.null(x[[1]]$delta))if(x[[1]]$delta!=0)paste("  delta: ",round(x[[1]]$delta,digits=2)),if(!is.null(x[[1]]$cp))if(x[[1]]$cp!=0)paste("  cp: ",round(x[[1]]$cp,digits=2))))
    
    cat("\n")
    mat <- matrix(c(round(x[[1]]$b,digits=3),round(x[[1]]$a,digits=3),round(x[[1]]$t,digits=3)),ncol=length(x[[1]]$t),byrow=TRUE)
    colnames(mat) <- rep("",length(x[[1]]$t))
    rownames(mat) <- c("Upper bounds","Lower bounds","Information fraction")
    print(mat)
    
    cat("\n");
    if(!is.null(x[[1]]$alab)) {
        mat_als <- matrix(c(round(x[[1]]$als,digits=3),round(x[[1]]$alab,digits=3)),ncol=length(x[[1]]$t),byrow=TRUE)
        colnames(mat_als) <- rep("",length(x[[1]]$t))
        rownames(mat_als) <- c("als","alab")
        print(mat_als)
        cat("\n")
    }

    cat("\n\n")
    cat(paste("interim data: \n"))
    cat(paste("\n\tT: ", x[[2]]$T, "\t z: ", round(x[[2]]$z, digits=3)))
    
    cat(paste("\n\nSecondary trial: \n\n"))
    cat(paste(length(x[[3]]$t)," stage group sequential design"))
    cat(paste("\n crp: ",round(x[[3]]$al,digits=3),"  SF: ",x[[3]]$SF,"  phi: ",x[[3]]$phi,if(!is.null(x[[3]]$Imax))paste("  Imax: ",round(x[[3]]$Imax,digits=2)),if(!is.null(x[[3]]$delta))if(x[[3]]$delta!=0)paste("  delta: ",round(x[[3]]$delta,digits=2)),if(!is.null(x[[3]]$cp))if(x[[3]]$cp!=0)paste("  cp: ",round(x[[3]]$cp,digits=2))))
    
    cat("\n")
    mat <- matrix(c(round(x[[3]]$b,digits=3),round(x[[3]]$a,digits=3),round(x[[3]]$t,digits=3)),ncol=length(x[[3]]$t),byrow=TRUE)
    colnames(mat) <- rep("",length(x[[3]]$t))
    rownames(mat) <- c("Upper bounds","Lower bounds","Information fraction")
    print(mat)
    
    cat(paste("\n\nSecondary trial outcome: \n"))
    cat(paste("\n\tT: ",x[[4]]$T,"  z: ",round(x[[4]]$z,digits=3),"\n\n"))    
}


#' Function \code{as.AGST} builds an adaptive group sequential trial object
#' 
#' @title as Adaptive Group Sequential Trial
#'
#' @param pT object of the \code{class} \code{GSTobj}; primary trial design
#' @param iD interim data; a list with the variables \code{T} and \code{z}; list(T = stage of interim analysis, z = interim z-statistic)
#' @param sT object of the \code{class} \code{GSTobj}; secondary trial design
#' @param sTo secondary trial outcome; a list with the variables \code{T} and \code{z}; list(T = stage where trial stops, z = z-statistic at stage where trial stops)
#' @return Returns a list containing the \code{pT}, \code{iD}, \code{sT} and \code{sTo} with \code{class}=\code{AGSTobj}
#' @author Niklas Hack \email{niklas.hack@@meduniwien.ac.at} and Werner Brannath \email{werner.brannath@@meduniwien.ac.at} 
#' @seealso \code{\link{AGSTobj}} 
#' @examples
#' pT=plan.GST(K=3,SF=4,phi=-4,alpha=0.05,delta=6,pow=0.9,compute.alab=TRUE,compute.als=TRUE)
#' iD=list(T=1, z=1.090728)
#' swImax=0.0625
#' I2min=3*swImax
#' I2max=3*swImax
#'
#' sT=adapt(pT=pT,iD=iD,SF=1,phi=0,cp=0.8,theta=5,I2min,I2max,swImax)
#' sTo=list(T=2, z=2.393)
#' AGST <- as.AGST(pT=pT,iD=iD,sT=sT,sTo=sTo)
#' 
#' @keywords methods
#' @export
as.AGST <- function(pT, iD, sT, sTo=NULL) {
    if(is.null(sT)) {
        warning("No design adaptation and hence no secondary trial; a GST object is created.")
        as.GST(GSD=pT, GSDo=iD)
    } else {
        if(is.null(sTo)) AGST <- list(pT=pT, iD=iD, sT=sT, sTo=list(T=0, z=0))
        else AGST <- list(pT=pT, iD=iD, sT=sT, sTo=sTo)
        class(AGST) <- "AGSTobj"
        AGST
    }
}
