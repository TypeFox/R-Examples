CO <- function(h, pT, x=pT$b[length(pT$t)], iD) {
    if(h == 0) k <- length(pT$t)
    else {
        if(is.null(pT$alab)) pT$alab <- comp.alab(GSD=pT)
        k <- j.alab(h, GSD=pT)
    }
    if(k-iD$T <= 0) {
        print("Cannot compute CO for T>=k")
        0
    } else {
        if(x == Inf) A(h, k, pT, iD)
        else {
            ub <- pbounds(h=h, pT=list(t=pT$t[1:k], b=c(pT$b[1:(k-1)], x), Imax=pT$Imax), iD=iD)
            if(k - iD$T == 1) 1-pnorm(ub)
            else {
                seqmon(a=pT$a[(iD$T+1):k],
                       b=ub[1:(k-iD$T)], t=pT$t[(iD$T+1):k]*pT$Imax-pT$t[iD$T]*pT$Imax,
                       int=500*array(c(1), k-iD$T))[2*(k-iD$T)]
            }
        }
    }
}


cerr <- function(h=0, pT, iD, bhh=NULL){
    if(is.null(pT$alab)) pT$alab <- comp.alab(GSD=pT)
    if(h == 0) CO(h=0, pT=pT, iD=iD)
    else {
        if(h %in% pT$alab[1:(length(pT$t)-1)]) A(h, k=j.alab(h, GSD=pT)+1, pT, iD)
        else {
            if(is.null(bhh)) bhh <- bh(h, pT);
            CO(h=h, pT=pT, x=bhh, iD=iD) }
    }
}

#' Calculates the conditional type I error rate of a GSD
#'
#' @title Conditional type I error rate (also called conditional rejection probability)
#' @param pT object of the \code{class} \code{GSTobj}; primary trial design
#' @param iD interim data; a list with the variables \code{T} and \code{z}; list(T = stage of interim analysis, z = interim z-statistic)
#' @return cer conditional type I error rate
#' @references Mueller, HH, Schaefer, H (2001) ''Adaptive group sequential design for clinical trials: Combining the advantages of adaptive and of classical group sequential approaches'', \emph{Biometrics}, 57, 886-891. 
#' @author Niklas Hack \email{niklas.hack@@meduniwien.ac.at} and Werner Brannath \email{werner.brannath@@meduniwien.ac.at}
#' @seealso \code{\link{plan.GST}}
#' @export
#' @examples
#' ##The following calculates the conditional type I error rate
#' ##under the null hypotesis after an adaptation at the second stage
#' ##of the primary trial.
#' pT=plan.GST(K=4,SF=1,phi=0,alpha=0.025,delta=6,pow=0.8,compute.alab=TRUE,compute.als=TRUE)
#' cer(pT=pT,iD=list(T=2, z=1.09))
#'
#' @keywords methods
cer <- function(pT, iD) {
    ##pT$t<-pT$t*pT$Imax
    cerr(h=0, pT, iD, bhh=NULL)
}


cer.r <- function(u, K, pT, iD) {
    if(K-iD$T <= 0) {
        print("Cannot compute CO for T>=K")
        0
    }
    b <- compBounds(t=1:K/K, t2 = pT$t*pT$Imax, iuse = pT$SF, asf = NULL,
                    alpha = u, phi = ifelse(is.null(pT$phi),0,pT$phi),
                    ztrun = 8)#$upper.bounds
    
    ub <- pbounds(h=0,pT=list(t=pT$t[1:K],b=b,Imax=pT$Imax),iD=iD)
    if(K-iD$T == 1) 1-pnorm(ub)
    else {
        seqmon(a=pT$a[(iD$T+1):K],
               b=ub[1:(K-iD$T)],t=pT$t[(iD$T+1):K]*pT$Imax-pT$t[iD$T]*pT$Imax,
               int=500*array(c(1),K-iD$T))[2*(K-iD$T)]
    }
}

rCER <- function(h=0, pT, iD, level=NULL) {
    pT$K <- length(pT$t)
    if(!is.null(level) && (pT$al != level)) {
        pT$al <- level
        pT$b  <- compBounds(t=1:pT$K/pT$K, t2 = pT$t*pT$Imax, iuse = pT$SF, asf = NULL, 
                            alpha = level, phi = ifelse(is.null(pT$phi),0,pT$phi), 
                            ztrun = 8)#$upper.bounds
    }
    CO(h=0, pT=pT, iD=list(T=iD$T, z=iD$z-h*sqrt(pT$t[iD$T]*pT$Imax)))
}
