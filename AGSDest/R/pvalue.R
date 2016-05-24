P.so.gsd <- function(h=0, GSD, GSDo) {
    if((GSDo$T < length(GSD$t)) && (GSDo$z < GSD$b[GSDo$T])) stop("stopping rule NOT met")
    else {
        if(GSDo$T == 1) 1-pnorm(GSDo$z-h*sqrt(GSD$t[1]*GSD$Imax))
        else {
            seqmon(a=GSD$a[1:GSDo$T],
                   b=c(GSD$b[1:(GSDo$T-1)], GSDo$z) - h*sqrt(GSD$t[1:GSDo$T]*GSD$Imax),
                   t=GSD$t[1:GSDo$T]*GSD$Imax, int=500*rep.int(1, GSDo$T))[2*GSDo$T]
        }
    }
}

P.so.ad <- function(h=0, pT, iD, sT, sTo, prec=0.0001) {
    cer <- sapply((iD$T+1):length(pT$a), function(j) {
        pT_new <- list(a=pT$a[1:j], b=pT$b[1:j], t=pT$t[1:j], Imax=pT$Imax)
        ifelse(j > iD$T, CO(h=0, pT=pT_new, iD=iD), 0)
    })
    
    P2 <- P.so.gsd(h=0, sT, sTo)
    k <- sum(ifelse(P2 > cer, 1, 0)) + iD$T + 1
    
    bkl <- qnorm(1 - pT$als[k])
    bku <- 1 - qnorm((pT$als[k] + pT$als[k-1])/2 + pT$als[k-1])
    
    pT_new <- list(a=pT$a[1:k], b=pT$b[1:k], t=pT$t[1:k], al=pT$als[k], Imax=pT$Imax)
    pT_new$b[k] <- bkl  
    ceu <- cerr(h=0, pT=pT_new, iD)
    u <- 0
    
    while(P2 > ceu && u == 0) {
        pT_new$b[k] <- pT_new$b[k] + (bkl+bku)/2
        bku <- bkl
        bkl <- bkl+(bkl+bku)/2
        ceu <- cerr(h=0, pT=pT_new, iD)
        if(abs(sword(h=0, k, pT, bku) - pT$als[k+1]) <= prec) {
            u <- sword(h=0, k, pT, bku)
            return(u)
        }
    }
    
    bku <- uniroot(f=function(bk) cerr(h=0, pT=list(a=pT$a[1:k], b=c(pT$b[1:k-1],bk), t=pT$t[1:k], al=pT$als[k], Imax=pT$Imax), iD)-P2, interval=c(bkl,bku))$root
    sword(h=0, k, pT, bku)
    
    
}

P.r.gsd <- function(h=0, GSD, GSDo, prec=0.001) {
    K <- length(GSD$t)
    GSD$K <- length(GSD$t)
    if(is.null(GSD$SF) || GSD$SF == 5) {
        b <- GSDo$z - h*sqrt(GSD$t[GSDo$T]*GSD$Imax)
        if(K == 1) 1-pnorm(b)
        else {
            seqmon(a=GSD$a,
                   b=GSD$b*(b/GSD$b[GSDo$T]),
                   t=(GSD$t*GSD$Imax)/(GSD$t[K]*GSD$Imax),int=500*rep.int(1,K))[2*K]
        }
    } else {
        if(GSD$SF == 7) {
            h1 <- (b - GSD$b[GSDo$T]) / sqrt(GSD$t[GSDo$T]*GSD$Imax)
            seqmon(a=GSD$a,
                   b=GSD$b+h1*sqrt(GSD$t),
                   t=(GSD$t*GSD$Imax)/(GSD$t[K]*GSD$Imax),int=500*rep.int(1, K))[2*K]
        } else {
            bisearch(function(u) GSDo$z-h*sqrt(GSD$t[GSDo$T]*GSD$Imax)-
                         compBounds(t=1:GSD$K/GSD$K, t2 = GSD$t*GSD$Imax, iuse = GSD$SF, asf = NULL, 
                                    alpha = u, phi = ifelse(is.null(GSD$phi),0,GSD$phi), 
                                    ztrun = 8)[GSDo$T],
                     c(0,1),signfl=-1, signfu=1,tol=prec)$root
        }
    }
}

P.r.ad <- function(h=0, pT, iD, sT, sTo, prec=0.001) {
    K <- length(pT$a)
    P2 <- P.r.gsd(h=0, GSD=sT, GSDo=sTo)
    bisearch(f=function(u) cer.r(u, K, pT, iD) - P2, c(0,1), signfl=-1, signfu=1,tol=prec)$root    
}

#' Calculates the repeated or stage-wise adjusted p-value of a GSD or a AGSD
#' 
#' @title Calculates the p-value
#' @param object object of the \code{class} \code{GSTobj} or of the \code{class} \code{AGSTobj}
#' @param type p-value type: repeated "r", stage-wise ordering "so" or both "b" (default: "b")
#' @details
#'  \code{object} can be an object of the \code{class} \code{GSTobj} or an object of the \code{class} \code{AGSTobj}. 
#'  The function identifies the \code{class} of the object and calculates the corresponding p-value (classical or adaptive).
#'  
#'  If \code{object} has \code{class} \code{GSTobj}, then a p-value for a classical GSD is calculated.
#'  \code{type} defines the type of confidence interval that is calculated
#'  \tabular{ll}{
#'    \code{"r"}\tab Repeated p-value for a classical GSD\cr
#'    \code{"so"}\tab Stage-wise adjusted p-value for a classical GSD\cr
#'  }
#'  If \code{object} has \code{class} \code{AGSTobj}, then a p-value for a GSD with design adaptation is calculated.
#'  \code{type} defines the type of confidence interval that is calculated 
#'  \tabular{ll}{
#'    \code{"r"}\tab Repeated p-value for a GSD with design adaptations\cr
#'    \code{"so"}\tab Stage-wise adjusted p-value for a GSD with design adaptations\cr
#'  }  
#' @return
#'  The function \code{pvalue} returns according to the \code{object} the classical or adaptive p-value for the final stage.
#'  If the parameter value has the \code{class} \code{GSTobj} the classical p-value is calculated. If the 
#'  parameter value has the \code{class} \code{AGSTobj} the adaptive p-value is calculated.
#'  
#'  The calculated p-values are saved as:
#'  
#'  \item{pvalue.r}{ repeated p-value}
#'  \item{pvalue.so}{ stage-wise adjusted p-value}
#' @note
#' The stage-wise adjusted p-value can only be calculated at the stage where the trial stops and is only valid if the stopping rule is met.
#'
#' The repeated p-value can be calculated at every stage of the trial and
#' not just at the stage where the trial stops and is also valid if the stopping rule is not met. 
#'
#' For calculating the sequential p-values at stage \code{T} the user has to specify the outcome \code{GSDo} in the object \code{GSTobj}
#' or \code{sTo} (secondary trial outcome) in the object \code{AGSTobj}. A trial outcome is a list of the form 
#' \code{list=(T=stage of interim analysis, z = interim z-statistic)}; see the example below.
#' @references
#' Brannath, W, Mehta, CR, Posch, M (2008) ''Exact confidence bounds following
#' adaptive group sequential tests'', \emph{Biometrics} accepted.
#'
#' Jennison, C, Turnbull, BW (1989) ''Repeated confidence intervals for group
#' sequential clinical trials'', \emph{Contr. Clin. Trials}, 5, 33-45.
#'
#' Mehta, CR, Bauer, P, Posch, M, Brannath, W (2007) ''Repeated confidence
#' intervals for adaptive group sequential trials'', \emph{Statistics in Medicine}, 26, 5422-5433.
#' 
#' Mueller, HH, Schaefer, H (2001) ''Adaptive group sequential design for clinical
#' trials: Combining the advantages of adaptive and of classical group sequential
#' approaches'', \emph{Biometrics}, 57, 886-891.
#'
#' Tsiatis,AA, Rosner,GL, Mehta,CR (1984) ''Exact confidence intervals
#' following a group sequential test'', \emph{Biometrics}, 40, 797-804.
#' @author Niklas Hack \email{niklas.hack@@meduniwien.ac.at} and Werner Brannath \email{werner.brannath@@meduniwien.ac.at} 
#' @seealso \code{\link{AGSTobj}}, \code{\link{GSTobj}} 
#' @examples
#' ##The following calculates the repeated p-value of a group sequential trial
#'
#' \dontrun{
#' GSD=plan.GST(K=4,SF=1,phi=0,alpha=0.025,delta=6,pow=0.8,compute.alab=TRUE,compute.als=TRUE)
#' 
#' GST<-as.GST(GSD=GSD,GSDo=list(T=2, z=3.1))
#'
#' pvalue(GST,type="r")
#'
#' ##The stage-wise adjusted p-value of a group sequential trial is calculated by
#' pvalue(GST,type="so")
#'
#' ##The repeated p-value at the earlier stage T=1 where the trial stopping rule is not met.
#' pvalue(as.GST(GSD,GSDo=list(T=1,z=0.7)),type="r")
#' 
#' ##If the stage-wise adjusted p-value is calculated at this stage, 
#' ##the function returns an error message
#'
#' pvalue(as.GST(GSD,GSDo=list(T=1,z=0.7)),type="so")
#'
#' ##The repeated and the stage-wise adjusted p-value of a 
#' ##group sequential trial after a design adaptation is calculated by
#' 
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
#' pvalue(AGST)
#'
#' ##The repeated p-value at the earlier stage T=2 where the stopping rule is not met.
#'
#' pvalue(as.AGST(pT,iD,sT,sTo=list(T=2,z=1.7)),type="r")
#'
#' ##If the stage-wise adjusted p-value is calculated at this stage, 
#' ##the function returns an error message
#'
#' pvalue(as.AGST(pT,iD,sT,sTo=list(T=2,z=1.7)),type="so")
#' }
#' @keywords methods
#' @export
pvalue <- function(object, type=c("r", "so")) {

    res <- list()
    
    if(class(object) == "GSTobj") {
        GSD <- object$GSD
        GSDo <- object$GSDo       
        
        if(!is.null(GSDo$z) || !is.null(GSDo$T)) {
            if("r" %in% type) res$pvalue.r <- P.r.gsd(h=0,GSD,GSDo)
            if("so" %in% type) {
                if(GSDo$z < GSD$b[GSDo$T] & GSDo$T < GSD$K) {
                    cat("pvalue.so : z < b[T]; Stopping rule NOT met.\n")
                    NULL
                }
                else res$pvalue.so <- P.so.gsd(h=0,GSD,GSDo)
            }
        } else {
            print("interim data missing")
            NULL
        }        
    } else if(class(object) == "AGSTobj") {
        if("r" %in% type) res$pvalue.r <- P.r.ad(h=0, object$pT, object$iD, object$sT, object$sTo)
        if("so" %in% type) {
            if(object$sTo$z < object$sT$b[object$sTo$T] && object$sTo$T < object$sT$K) {
                cat("pvalue.so : z < b[T]; Stopping rule NOT met.\n")
                NULL
            }
            else res$pvalue.so <- P.so.ad(h=0, object$pT, object$iD, object$sT, object$sTo)
        }             
    }

    res
}
    
