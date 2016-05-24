#' The package allows to compute repeated confidence intervals as well as 
#' confidence intervals based on the stage-wise ordering in groug sequential designs (GSD; see Jennison and Turnbull, 1989; Tsiatis, Rosner, Mehta, 1984) and
#' adaptive groug sequential designs (Mehta, Bauer, Posch, Brannath, 2007; Brannath, Mehta, Posch, 2008). 
#' For adaptive group sequential designs the confidence intervals are based on the conditional 
#' rejection probability principle of Mueller and Schaefer (2001). This principle allows us 
#' to perform data dependent changes to the sample size, the spending function, 
#' and the number and spacing of interim looks while preserving the overall type I error rate. 
#' Currently the procedures do not support the use of futility boundaries as well as more than one adaptive interim analysis. 
#' Furthermore, the package is currently restricted to the computation of lower one-sided confidence intervals.
#'
#' Package: AGSDest
#'
#' Type: Package
#' Version: 2.2
#' Date: 2015-01-12
#' License: GPL Version 2 or later
#' 
#' Main functions:
#' \code{adapt}: Performs adaptations at an interim analysis of a GSD to the sample size, number of interim stages and spending function based on the conditional power in a GSD at an interim analysis; the result is a secondary trial
#' \code{plan.GST}: Plans a group sequential trial
#' \code{cer}: Computes the conditional type I error rate (also called conditional rejection probability) of a GSD at an interim analysis
#' \code{typeIerr}: Computes the type I error rate of a GSD
#' \code{pvalue}: Computes the repeated or stage-wise adjusted p-value for a classical GSD or for a GSD with design adaptations
#' \code{seqconfint}: Computes the repeated confidence bound and confidence bound based on the stage-wise ordering for a GSD or for a GSD with design adaptations
#'
#' Subfunctions:
#' \code{as.GST}: Builds a group sequential trial object
#' \code{as.AGST}: Builds an adaptive group sequential trial object
#'
#' @docType package
#' @title Estimation in adaptive group sequential trials
#' @author Niklas Hack \email{niklas.hack@@meduniwien.ac.at} and Werner Brannath \email{werner.brannath@@meduniwien.ac.at}
#' @references
#' Brannath, W, Mehta, CR, Posch, M (2008) ''Exact confidence bounds following adaptive group sequential tests'', \emph{Biometrics} accepted.
#' 
#' Jennison, C, Turnbull, BW (1989) ''Repeated confidence intervals for group sequential clinical trials'', \emph{Contr. Clin. Trials}, 5, 33-45.
#' 
#' Mehta, CR, Bauer, P, Posch, M, Brannath, W (2007) ''Repeated confidence intervals for adaptive group sequential trials'', \emph{Statistics in Medicine}, 26, 5422-5433.
#' 
#' Mueller, HH, Schaefer, H (2001) ''Adaptive group sequential design for clinical trials: Combining the advantages of adaptive and of classical group sequential approaches'', \emph{Biometrics}, 57, 886-891.
#' 
#' O'Brien, PC, Fleming, TR (1979) ''A multiple testing procedure for clinical trials'', \emph{Biometrics}, 35 , 549-556
#' 
#' Schoenfeld, D (2001) ''A simple Algorithm for Designing Group Sequential Clinical Trials'', \emph{Biometrics}, 27, 972-974
#' 
#' Tsiatis,AA, Rosner,GL, Mehta,CR (1984) ''Exact confidence intervals following a group sequential test'', \emph{Biometrics}, 40, 797-804.
#' 
#' @examples
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
#' sTo=list(T=2, z=2.393)
#' AGST <- as.AGST(pT=pT,iD=iD,sT=sT,sTo=sTo)
#' 
#' ##The following calculates the stage-wise adjusted p-value
#' ##of a group sequential trial after a design adaptation
#' pvalue(AGST,type="so")
#' 
#' ##and the corresponding confidence bound based on the stage-wise ordering.
#' seqconfint(AGST,type="so")
#' 
#' ##Both, the p-value and the confidence interval can be calculated by
#' ##the summary function
#' \dontrun{
#' summary(AGST,ctype="so",ptype="so")
#' }
#' @keywords methods list datasets
#'
#' @import stats graphics grDevices ldbounds
#' @useDynLib AGSDest mainf
#' @name AGSDest
NULL
