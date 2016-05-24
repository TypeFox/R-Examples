#' The \code{GSTobj} includes design and outcome of primary trial.
#'
#' A \code{GSTobj} object is designed.
#'
#' The function \code{summary} returns an object of \code{class} \code{GSTobj}.
#'
#' \code{ctype} defines the type of confidence interval that is calculated.
#' \tabular{ll}{
#'    \code{"r"}\tab Repeated confidence bound for a classical GSD\cr
#'    \code{"so"}\tab Confidence bound for a classical GSD based on the stage-wise ordering\cr
#'  }
#'  The calculated confidence bounds are saved as:
#'  \tabular{ll}{
#'  \code{cb.r}\tab repeated confidence bound\cr
#'  \code{cb.so}\tab confidence bound based on the stage-wise ordering\cr
#'  }
  
#'  \code{ptype} defines the type of p-value that is calculated.
#'    \tabular{ll}{
#'    \code{"r"}\tab Repeated p-value for a classical GSD\cr
#'    \code{"so"}\tab Stage-wise adjusted p-value for a classical GSD \cr
#'  }
  
#'  The calculated p-values are saved as:
#'  \tabular{ll}{
#'  \code{pvalue.r}\tab repeated p-value\cr
#'  \code{pvalue.so}\tab stage-wise adjusted p-value\cr
#'  }
  
#'  \code{etype} defines the type of point estimate
#'    \tabular{ll}{
#'    \code{"ml"}\tab maximum likelihood estimate (ignoring the sequential nature of the design)\cr
#'    \code{"mu"}\tab median unbiased estimate (stage-wise lower confidence bound at level 0.5) for a classical GSD\cr
#'    \code{"cons"}\tab Conservative estimate (repeated lower confidence bound at level 0.5) for a classical GSD\cr   
#'    }
  
#'  The calculated point estimates are saved as:
#'  \tabular{ll}{
#'  \code{est.ml}\tab Maximum likelihood estimate\cr
#'  \code{est.mu}\tab Median unbiased estimate\cr
#'  \code{est.cons}\tab Conservative estimate\cr
#'  }
  
#'  The stage-wise adjusted confidence interval and p-value and the median unbiased point estimate can only be calculated at the stage where the trial stops and is only valid if the stopping rule is met.
#'
#' The repeated confidence interval and repeated p-value, conservative estimate and maximum likelihood estimate can be calculated at every stage of the trial and
#' not just at the stage where the trial stops and is also valid if the stopping rule is not met.
#' For calculating the repeated confidence interval or p-value at any stage of the trial the user has to specify the outcome \code{GSDo} in the object \code{GSTobj} (see example below).
#' 
#' @title Group sequential trial object (GSTobj)
#'
#' @param x object of the \code{class} \code{GSTobj}
#' @param object object of the \code{class} \code{GSTobj}
#' @param main Title of the plots (default: "GSD")
#' @param print.pdf option; if TRUE a pdf file is created. Instead of setting print.pdf to TRUE, the user can specify a character string giving the name or the path of the file.
#' @param ctype confidence type: repeated "r" or stage-wise ordering "so" (default: c("r", "so"))
#' @param ptype p-value type: repeated "r" or stage-wise ordering "so" (default: c("r", "so"))
#' @param etype point estimate: maximum likelihood "ml", median unbiased "mu" or conservative "cons" (default: c("ml", "mu", "cons"))
#' @param overwrite option; if TRUE all old values are deleted and new values are calculated (default: FALSE)
#' @param ... additional arguments.
#' 
#' @return An object of class \code{GSTobj}, is basically a list with the elements
#' \item{cb.so}{ confidence bound based on the stage-wise ordering}
#' \item{cb.r}{ repeated confidence bound}
#' \item{pvalue.so}{ stage-wise adjusted p-value}
#' \item{pvalue.r}{ repeated p-value}
#' \item{est.ml}{ maximum likelihood estimate}
#' \item{est.mu}{ median unbiased point estimate}
#' \item{est.cons}{ conservative point estimate}
#' \item{GSD}{}
#' \item{K}{ number of stages}
#' \item{al}{ alpha (type I error rate)}
#' \item{a}{ lower critical bounds of group sequential design (are currently always set to -8)}
#' \item{b}{ upper critical bounds of group sequential design }
#' \item{t}{ vector with cumulative information fraction}
#' \item{SF}{ spending function (for details see below)}
#' \item{phi}{ parameter of spending function when SF=3 or 4 (for details see below)}
#' \item{alab}{ alpha-absorbing parameter values of group sequential design}
#' \item{als}{ alpha-values ''spent'' at each stage of group sequential design}
#' \item{Imax}{ maximum information number}
#' \item{delta}{ effect size used for planning the primary trial}
#' \item{cp}{ conditionla power of the trial}
#' \item{GSDo}{}
#' \item{T}{ stage where trial stops}
#' \item{z}{ z-statistic at stage where trial stops}
#' @author Niklas Hack \email{niklas.hack@@meduniwien.ac.at} and Werner Brannath \email{werner.brannath@@meduniwien.ac.at}
#' @note
#' \code{SF} defines the spending function.
#' \code{SF} = 1 O'Brien and Fleming type spending function of Lan and DeMets (1983)
#' \code{SF} = 2 Pocock type spending function of Lan and DeMets (1983)
#' \code{SF} = 3 Power family (\eqn{c_\alpha* t^\phi}). phi must be greater than 0.
#' \code{SF} = 4 Hwang-Shih-DeCani family.\eqn{(1-e^{-\phi t})/(1-e^{-\phi})}, where \code{phi} cannot be 0.
#' A value of \code{SF}=3 corresponds to the power family. Here, the spending function is \eqn{t^{\phi}},
#' where phi must be greater than 0. A value of \code{SF}=4 corresponds to the Hwang-Shih-DeCani family,
#' with the spending function \eqn{(1-e^{-\phi t})/(1-e^{-\phi})}, where phi cannot be 0.
#' If a path is specified for \code{print.pdf}, all \ must be changed to /. If a filename is specified the ending of the file must be (.pdf).
#' In the current version \code{a} should be set to rep(-8,K)
#'
#' @seealso \code{\link{GSTobj}}, \code{\link{print.GSTobj}}, \code{\link{plot.GSTobj}}, \code{\link{summary.GSTobj}}
#' @examples
#' GSD=plan.GST(K=4,SF=1,phi=0,alpha=0.025,delta=6,pow=0.8,compute.alab=TRUE,compute.als=TRUE)
#' GST<-as.GST(GSD=GSD,GSDo=list(T=2, z=3.1))
#' GST
#' plot(GST)
#' GST<-summary(GST)
#' plot(GST)
#' ##The repeated confidence interval, p-value and maximum likelihood estimate
#' ##at the earlier stage T=1 where the trial stopping rule is not met.
#' summary(as.GST(GSD,GSDo=list(T=1,z=0.7)),ctype="r",ptype="r",etype="ml")
#' \dontrun{
#' ##If e.g. the stage-wise adjusted confidence interval is calculated at this stage,
#' ##the function returns an error message
#' summary(as.GST(GSD,GSDo=list(T=1,z=0.7)),ctype="so",etype="mu")
#' }
#' @keywords datasets
#' @export
GSTobj <- function(x, ...) UseMethod("GSTobj")


#' @rdname GSTobj
#' @export
plot.GSTobj <- function(x, main="GSD", print.pdf=FALSE, ...) {

                                        #par.old=par(no.readonly = TRUE)
                                        #par(mfrow=c(1,1),oma = c(0,12,0,0))

    if(is.null(x)) {
        cat(paste("No design adaptation performed.\n"))
        return(0)
    } else{
        if(!is.null(x$t)) {
            
            if(mode(print.pdf) == "logical" && print.pdf == FALSE) {
                nf <- layout(matrix(c(1,2),1,2,byrow=TRUE), c(2.2,3), c(3,3))
                layout.show(nf)
            }
            if(mode(print.pdf) == "logical" && print.pdf == TRUE) pdf("GST.pdf")
            if(mode(print.pdf) == "character") pdf(print.pdf)
            
            plot(1, 1, type="n", ann=FALSE, axes=FALSE, col="black")
            box(col="black")
            
            text(1,1,paste(length(x$t)," stage GSD,\n",switch(x$SF,"OBF","Pocock","Power family","Hwang-Shih-DeCani","self specified","NA")," boundaries \nat level ",round(x$al,digits=3),
                           if(!is.null(x$Imax))paste(",\nImax=",round(x$Imax,digits=2)),
                           if(!is.null(x$delta))if(x$delta!=0)paste(",\ndelta=",round(x$delta,digits=2)),
                           paste("\n"),
                           if(!is.null(x$pvalue.r))paste("\npvalue.r = ",round(x$pvalue.r,digits=3)),
                           if(!is.null(x$pvalue.so))paste("\npvalue.so = ",round(x$pvalue.so,digits=3)),
                           if(!is.null(x$cb.r))paste("\ncb.r = ",round(x$cb.r,digits=3)),
                           if(!is.null(x$cb.so))paste("\ncb.so = ",round(x$cb.so,digits=3)),
                           if(!is.null(x$est.ml))paste("\nest.ml = ",round(x$est.ml,digits=3)),
                           if(!is.null(x$est.mu))paste("\nest.mu = ",round(x$est.mu,digits=3))),cex=0.8)

            
            plot(x$t,x$b,axes=FALSE,ylim=c(0,max(x$b)+1.5),xlab="Cumulative Information Fraction",ylab="Wald Teststatistic",cex.lab=1.25,type="l",
                 main=main, lwd=1.5,cex.main=0.9,cex.lab=0.9,cex=0.9,cex.axis=0.9,cex.sub=0.9,...)
            
            axis(2,cex.axis=0.9)
            axis(1, at=x$t, labels=(round(x$t,digits=3)),cex.axis=0.9)
            points(x$t,x$b,lt=1,pch=19)
            
                                        #legend(x$t[1],min(x$b),c("Boundaries"),lty=c(1),pch=c(19),cex=0.8,...)
            

            if(mode(print.pdf) == "logical" && print.pdf == TRUE) {
                dev.off()
                cat("File GST.pdf is created\n")
            }
            if(mode(print.pdf) == "character") {
                dev.off()
                cat("File", print.pdf, "is created\n")
            }
        }
    }

    if(is.null(x[[1]])) cat(paste("No design adaptation performed.\n"))
    else {
        if(!is.numeric(x[[1]])) {
            if(!is.null(x[[1]]$t)) {
                
                if(mode(print.pdf) == "logical" && print.pdf == FALSE) {
                    nf <- layout(matrix(c(1,2), 1, 2, byrow=TRUE), c(2.2,3), c(3,3))
                    layout.show(nf)
                }
                if(mode(print.pdf) == "logical" && print.pdf == TRUE) pdf("GST.pdf")
                if(mode(print.pdf) == "character") pdf(print.pdf)
                
                plot(1, 1, type="n", ann=FALSE, axes=FALSE, col="black")
                box(col="black")

                text(1,1,paste(length(x[[1]]$t)," stage GSD,\n",switch(x[[1]]$SF,"OBF","Pocock","Power family","Hwang-Shih-DeCani","self specified","NA")," boundaries \nat level ",round(x[[1]]$al,digits=3),
                               if(!is.null(x[[1]]$Imax))paste(",\nImax=",round(x[[1]]$Imax,digits=2)),
                               if(!is.null(x[[1]]$delta))if(x[[1]]$delta!=0)paste(",\ndelta=",round(x[[1]]$delta,digits=2)),
                               if(!is.null(x[[2]]$T))paste("\n\nTrial stops at look T = ",x[[2]]$T,"\nwith z = ",round(x[[2]]$z,digits=3)),
                               paste("\n"),
                               if(!is.null(x$pvalue.r))paste("\npvalue.r = ",round(x$pvalue.r,digits=3)),
                               if(!is.null(x$pvalue.so))paste("\npvalue.so = ",round(x$pvalue.so,digits=3)),
                               if(!is.null(x$cb.r))paste("\ncb.r = ",round(x$cb.r,digits=3)),
                               if(!is.null(x$cb.so))paste("\ncb.so = ",round(x$cb.so,digits=3)),
                               if(!is.null(x$est.ml))paste("\nest.ml = ",round(x$est.ml,digits=3)),
                               if(!is.null(x$est.mu))paste("\nest.mu = ",round(x$est.mu,digits=3))),cex=0.8)

                
                plot(x[[1]]$t,x[[1]]$b,axes=FALSE,ylim=c(0,max(x[[1]]$b)+1.5),xlab="Cumulative Information Fraction",ylab="Wald Teststatistic",cex.lab=1.25,type="l",
                     main=main, lwd=1.5,cex.main=0.9,cex.lab=0.9,cex=0.9,cex.axis=0.9,cex.sub=0.9,...)
                
                axis(2,cex.axis=0.9)
                axis(1, at=x[[1]]$t, labels=(round(x[[1]]$t,digits=3)),cex.axis=0.9)
                points(x[[1]]$t,x[[1]]$b,lt=1,pch=19)

                if(!is.null(x[[2]]$T)||x[[2]]$T == 0) {
                    points(x[[1]]$t[x[[2]]$T],x[[2]]$z,pch=21)
                                        #legend(x[[1]]$t[1],min(x[[1]]$b,x[[2]]$z)-0.3,c("Standardized test statistic","Boundaries"),lty=c(0,1),pch=c(21,19),cex=0.8,...)
                }
                                        #else legend(x[[1]]$t[1],min(x[[1]]$b),c("Boundaries"),lty=c(1),pch=c(19),cex=0.8,...)


                if(mode(print.pdf) == "logical" && print.pdf == TRUE) {
                    dev.off()
                    cat("File GST.pdf is created\n")
                }
                if(mode(print.pdf) == "character") {
                    dev.off()
                    cat("File", print.pdf, "is created\n")
                }
            }
        }
    }
}

#' @rdname GSTobj
#' @export
print.GSTobj <- function(x, ...) {

    if(!is.null(x$cb.r)) cat(paste("repeated lower confidence bound: ",round(x$cb.r,digits=3),"\n\n"))
    if(!is.null(x$cb.so)) cat(paste("stage-wise adjusted lower confidence bound: ",round(x$cb.so,digits=3),"\n\n"))
    if(!is.null(x$pvalue.r)) cat(paste("repeated p-value: ",round(x$pvalue.r,digits=3),"\n\n"))
    if(!is.null(x$pvalue.so)) cat(paste("stage-wise adjusted p-value: ",round(x$pvalue.so,digits=3),"\n\n"))

    if(!is.null(x$est.ml)) cat(paste("maximum likelihood estimate: ",round(x$est.ml,digits=3),"\n\n"))
    if(!is.null(x$est.mu)) cat(paste("median unbiased estimate: ",round(x$est.mu,digits=3),"\n\n"))
    if(!is.null(x$est.cons)) cat(paste("conservative estimate: ",round(x$est.cons,digits=3),"\n\n"))

    if(is.null(x)) cat(paste("No design adaptation performed.\n"))
    else{
        if(!is.null(x$t)) {

            cat(paste(length(x$t)," stage group sequential design"))
            cat(paste("\n",expression(alpha),": ",round(x$al,digits=3),"  SF: ",x$SF,"  phi: ",x$phi,if(!is.null(x$Imax))paste("  Imax: ",round(x$Imax,digits=2)),if(!is.null(x$delta))if(x$delta!=0)paste("  delta: ",round(x$delta,digits=2)),if(!is.null(x$cp))if(x$cp!=0)paste("  cp: ",round(x$cp,digits=2))))

            cat("\n")
            mat=matrix(c(round(x$b,digits=3),round(x$a,digits=3),round(x$t,digits=3)),ncol=length(x$t),byrow=TRUE)
            colnames(mat)=rep("",length(x$t))
            rownames(mat)=c("Upper bounds","Lower bounds","Information fraction")
            print(mat)

            cat("\n");
            if(!is.null(x$alab)){
                mat_als=matrix(c(round(x$als,digits=3),round(x$alab[1:length(x$t)],digits=3)),ncol=length(x$t),byrow=TRUE)
                colnames(mat_als)=rep("",length(x$t))
                rownames(mat_als)=c("als","alab")
                print(mat_als)
                cat("\n\n")
            }
        }
    }


    if(is.null(x[[1]])) cat(paste("No design adaptation performed.\n"))
    else {

        if(!is.numeric(x[[1]])) {
            if(!is.null(x[[1]]$t)) {

                cat(paste(length(x[[1]]$t)," stage group sequential design"))
                cat(paste("\n",expression(alpha),": ",round(x[[1]]$al,digits=3),"  SF: ",x[[1]]$SF,"  phi: ",x[[1]]$phi,if(!is.null(x[[1]]$Imax))paste("  Imax: ",round(x[[1]]$Imax,digits=2)),if(!is.null(x[[1]]$delta))if(x[[1]]$delta!=0)paste("  delta: ",round(x[[1]]$delta,digits=2)),if(!is.null(x[[1]]$cp))if(x[[1]]$cp!=0)paste("  cp: ",round(x[[1]]$cp,digits=2))))

                cat("\n")
                mat=matrix(c(round(x[[1]]$b,digits=3),round(x[[1]]$a,digits=3),round(x[[1]]$t,digits=3)),ncol=length(x[[1]]$t),byrow=TRUE)
                colnames(mat)=rep("",length(x[[1]]$t))
                rownames(mat)=c("Upper bounds","Lower bounds","Information fraction")
                print(mat)

                if(!is.null(x[[1]]$alab)) {
                    mat_als=matrix(c(round(x[[1]]$als,digits=3),round(x[[1]]$alab[1:length(x[[1]]$t)],digits=3)),ncol=length(x[[1]]$t),byrow=TRUE)
                    colnames(mat_als)=rep("",length(x[[1]]$t))
                    rownames(mat_als)=c("als","alab")
                    print(mat_als)
                    cat("\n")
                }

                if(!is.null(x[[2]]$T)) {
                    cat("\ngroup sequential design outcome:\n")
                    cat(paste("\n\tT: ",x[[2]]$T,"  z: ",round(x[[2]]$z,digits=3),"\n\n"))
                }
            }
        }
    }
}

#' @rdname GSTobj
#' @export
summary.GSTobj <- function(object, ctype=c("r", "so"), ptype=c("r", "so"), etype=c("ml", "mu", "cons"), overwrite=FALSE,...) {

    GSD <- object$GSD
    GSDo <- object$GSDo
    
    if(!is.null(GSDo$z) && !is.null(GSDo$T)) {
        if("r" %in% ctype) {
            if(is.null(object$cb.r) || object$cb.r==0 || (!is.null(object$cb.r) && overwrite)) {
                object$cb.r <- cb.r.gsd(GSD, GSDo)
            }
        }
        
        if("so" %in% ctype) {
            if(is.null(object$cb.so) || object$cb.so==0 || (!is.null(object$cb.so) && overwrite)) {
                if(GSDo$z < GSD$b[GSDo$T] && GSDo$T < GSD$K) cat("cb.so : z < b[T]; Stopping rule NOT met.\n")
                else object$cb.so <- cb.so.gsd(GSD, GSDo)
            }
        }

        if("r" %in% ptype) {
            if(is.null(object$pvalue.r) || object$pvalue.r==0 || (!is.null(object$pvalue.r) && overwrite)) {
                object$pvalue.r=P.r.gsd(h=0,GSD,GSDo)
            }
        }
        
        if("so" %in% ptype) {
            if(is.null(object$pvalue.so) || object$pvalue.so==0 || (!is.null(object$pvalue.so) && overwrite)) {
                if(GSDo$z < GSD$b[GSDo$T] && GSDo$T < GSD$K) {
                    cat("pvalue.so : z < b[T]; Stopping rule NOT met.\n")
                }
                else object$pvalue.so <- P.so.gsd(h=0, GSD, GSDo)
            }
        }
        
        if("ml" %in% etype) {
            if(is.null(object$est.ml) || object$est.ml==0 || (!is.null(object$est.ml) && overwrite)) {
                object$est.ml <- GSDo$z / sqrt(GSD$t[GSDo$T]*GSD$Imax)
            }
        }
        
        if("mu" %in% etype) {
            if(is.null(object$est.mu) || object$est.mu==0 || (!is.null(object$est.mu) && overwrite)){
                if(GSDo$z < GSD$b[GSDo$T]) cat("est.mu : z < b[T]; Stopping rule NOT met.\n")
                else object$est.mu <- cb.so.gsd(GSD, GSDo, level=0.5)
            }
        }
        
        if("cons" %in% etype) {
            if(is.null(object$est.cons) || object$est.cons==0 || (!is.null(object$est.cons) && overwrite))
                object$est.cons <- cb.r.gsd(GSD, GSDo, level=0.5)
        }

        object
    } else print("Missing interim data")
}

#' @rdname GSTobj
#' @export
print.summary.GSTobj <- function(x, ...) {

    if(!is.null(x$type))cat(paste("type: ",x$type,"\n\n"))

    if(!is.null(x$cb.r)) cat(paste("repeated lower confidence bound: ",round(x$cb.r,digits=3),"\n\n"))
    if(!is.null(x$cb.so)) cat(paste("stage-wise adjusted lower confidence bound: ",round(x$cb.so,digits=3),"\n\n"))
    if(!is.null(x$pvalue.r)) cat(paste("repeated p-value: ",round(x$pvalue.r,digits=3),"\n\n"))
    if(!is.null(x$pvalue.so)) cat(paste("stage-wise adjusted p-value: ",round(x$pvalue.so,digits=3),"\n\n"))

    if(!is.null(x$est.ml)) cat(paste("maximum likelihood estimate: ",round(x$est.ml,digits=3),"\n\n"))
    if(!is.null(x$est.mu)) cat(paste("median unbiased estimate: ",round(x$est.mu,digits=3),"\n\n"))
    if(!is.null(x$est.cons)) cat(paste("conservative estimate: ",round(x$est.cons,digits=3),"\n\n"))

    if(!is.null(x$GSD)) {
        cat(paste(length(x$GSD$a)," stage group sequential design"))
        cat(paste("\n",expression(alpha),": ",round(x$GSD$al,digits=3),"  SF: ",x$GSD$SF,"  phi: ",x$GSD$phi,if(!is.null(x$Imax))paste("  Imax: ",round(x$Imax,digits=2)),if(!is.null(x$delta))if(x$delta!=0)paste("  delta: ",round(x$delta,digits=2)),if(!is.null(x$cp))if(x$cp!=0)paste("  cp: ",round(x$cp,digits=2))))

        cat("\n")
        mat=matrix(c(round(x$GSD$b,digits=3),round(x$GSD$a,digits=3),round(x$GSD$t,digits=3)),ncol=length(x$GSD$t),byrow=TRUE)
        colnames(mat)=rep("",length(x$GSD$t))
        rownames(mat)=c("Upper bounds","Lower bounds","Information fraction")
        print(mat)

        cat("\n");
        if(!is.null(x$GSD$alab)){
            mat_als=matrix(c(round(x$GSD$als,digits=3),round(x$GSD$alab,digits=3)),ncol=length(x$GSD$t),byrow=TRUE)
            colnames(mat_als)=rep("",length(x$GSD$t))
            rownames(mat_als)=c("als","alab")
            print(mat_als)
            cat("\n\n")
        }

        cat("\n\ngroup sequential design outcome:\n")
        cat(paste("\n\tT: ",x$GSDo$T,"  z: ",round(x$GSDo$z,digits=3),"\n\n"))
    }
}

#'
#' Plans a group sequential trial (\code{GST})
#' @title Plans a group sequential trial (GST)
#' @param K number of stages 
#' @param t vector with the cumulative information fraction (default: (1:K)/K) 
#' @param Imax maximum information number (default: NULL)
#' @param SF spending function (for details see below)
#' @param phi parameter of spending function when SF=3 or 4 (See below) 
#' @param alpha alpha (type I error rate) 
#' @param delta effect size (alternative)(default: NULL)
#' @param pow power (default: NULL)
#' @param compute.alab specify if alpha-absorbing parameter values should be calculated (default: TRUE) 
#' @param compute.als specify if alpha-values ''spent'' at every stage should be calculated (default: TRUE) 
#' @details
#' The user has to specify either \code{Imax} or \code{delta} and \code{pow}.
#' If all three items are specified, the pre-defined maximum information number is newly calculated from the information for \code{delta} and \code{power}, and \code{Imax} is overwritten.
#'
#' \code{SF} defines the spending function.
#' \tabular{ll}{
#' \code{SF} = \tab 1 O'Brien and Fleming type spending function of Lan and DeMets (1983)\cr
#' \code{SF} = \tab 2 Pocock type spending function of Lan and DeMets (1983)\cr
#' \code{SF} = \tab 3 Power family (\eqn{c_\alpha* t^\phi}); \code{phi} must be greater than 0\cr
#' \code{SF} = \tab 4 Hwang-Shih-DeCani family; \eqn{(1-e^{-\phi t})/(1-e^{-\phi})}, where \code{phi} cannot be 0\cr
#' }
#' @return
#' \code{plan.GST} returns an object of the \code{class} \code{GSTobj}.
#' An object of class \code{GSTobj} is a list containing the following
#' components:
#' \item{K}{ number of stages}
#' \item{a}{ lower critical bounds of group sequential design (are currently always set to -8)}
#' \item{b}{ upper critical bounds of group sequential design}
#' \item{t}{ vector with cumulative information fraction}
#' \item{al}{ alpha (type I error) }
#' \item{SF}{ spending function }
#' \item{phi}{ parameter of spending function when SF=3 or 4 (See below)}
#' \item{Imax}{ maximum information number}
#' \item{delta}{ effect size used for planning the primary trial}
#' @references
#' Brannath, W, Mehta, CR, Posch, M (2008) ''Exact confidence bounds following
#' adaptive group sequential tests'', \emph{Biometrics} accepted.
#'
#' @author Niklas Hack \email{niklas.hack@@meduniwien.ac.at} and Werner Brannath \email{werner.brannath@@meduniwien.ac.at}
#' @seealso \code{\link{GSTobj}}, \code{\link{print.GSTobj}}, \code{\link{plot.GSTobj}}
#' @examples
#' ##The following plans an O'Brien and Flaming group sequential design (GSD) 
#' ##with 4 stages and equally spaced looks.
#' pT <- plan.GST(K=4, SF=1, phi=0, alpha=0.025, delta=6, pow=0.8, compute.alab=TRUE, compute.als=TRUE)
#' 
#' @keywords methods
#' @export
plan.GST <- function(K, t=(1:K)/K, Imax=NULL, SF, phi, alpha, delta=NULL, pow=NULL, compute.alab=TRUE, compute.als=TRUE) {
    
    if(!is.null(pow) && (pow > 1)) pow <- pow/100

    if(SF==4 & phi == 0) stop("phi must be unequal 0")
    else {
        if((!is.null(delta) & !is.null(pow)) | !is.null(Imax)) {
            if(!is.null(delta) & !is.null(pow) & !is.null(Imax)) print("\nGST is planned only on the basis of Imax.\n")
            t <- t[order(t)]/max(t)
            if(is.null(Imax)) Imax <- 1
            b <- compBounds(t=t, t2 = t*Imax, iuse = SF, asf = NULL,
                            alpha = alpha, phi = phi,ztrun = 8)#$upper.bounds
            a <- rep(-8,K)
            
            if(!is.null(delta)) {  
                hmin <- qnorm(1-alpha)-qnorm(1-pow)
                hmax <- b[K]-qnorm(1-pow)

                if(K == 1) h <- hmax
                else {
                    h <- uniroot(function(x) seqmon(a=a[1:K],b=pbounds(h=x,pT=list(t=t,b=b,Imax=Imax),iD=list(T=0)),t=t[1:K],
                        int=500*array(c(1),K))[2*K]-pow,c(hmin,hmax))$root          
                }
                Imax <- (h/delta)^2
            }
            else delta <- 0
            b <- compBounds(t=(1:K)/K, t2 = Imax*t, iuse = SF, asf = NULL,
                            alpha = alpha, phi = phi,ztrun = 8)

            GSD <- list(K=K, al=alpha, t=t, SF=SF, phi=phi, a=a, b=b, Imax=Imax, delta=delta)
            
            if(compute.alab) GSD$alab <- comp.alab(GSD=GSD)          
            if(compute.als) GSD$als <- comp.als(GSD=GSD)           

            GSD$cp <- cp(GSD)
            class(GSD) <- "GSTobj"
            GSD
        }
        else stop("Please specify either delta and pow or Imax")
    }
}

#'
#' Function \code{as.GST} builds a group sequential trial object
#' 
#' @title as Group Sequential Trial
#'
#' @param GSD object of the \code{class} \code{GSTobj}; group sequential design
#' @param GSDo group sequential design outcome;  a list with the variables \code{T} and \code{z}; list(T = stage where trial stops, z = z-statistic at stage where trial stops)
#' @return Returns a list containing the \code{GSD} and \code{GSDo} with \code{class}=\code{GSTobj}
#' @author Niklas Hack \email{niklas.hack@@meduniwien.ac.at} and Werner Brannath \email{werner.brannath@@meduniwien.ac.at}
#' @seealso \code{\link{GSTobj}}
#' @examples 
#' GSD <- plan.GST(K=4,SF=1,phi=0,alpha=0.025,delta=6,pow=0.8,compute.alab=TRUE,compute.als=TRUE)
#' GSDo <- list(T=2, z=3.1)
#' GST <- as.GST(GSD=GSD,GSDo=GSDo)
#' GST
#' 
#' @keywords methods
#' @export
as.GST <- function(GSD, GSDo) {
    ##if(GSDo$z<GSD$b[GSDo[1][[1]]])print("Stopping rule not met")
    ##else{
    GST <- list(GSD=GSD, GSDo=GSDo)
    class(GST) <- "GSTobj"
    GST
}
