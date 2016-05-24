#' testGyriq
#'
#' Calculates the p-value of the kinship-adjusted SNP-set association test for
#' censored traits
#'
#' If the lower and upper bounds of sliding windows are not provided, the test 
#' is performed once on the whole SNP-set \code{G}. Otherwise, the score
#' statistic and the p-value are computed for each window sequentially.
#' 
#' In each run, the score statistic, which has a quadratic form following a 
#' mixture of chi-squared variables, is calculated from the completed vector of
#' residuals and a kernel matrix. The p-value is obtained using a permutation
#' approach based on matching moments described in Lee et al. (2012), a standard
#' permutation procedure or the Davies approximation (Davies, 1980) implemented 
#' in the package \pkg{CompQuadForm} (Duchesne and Lafaye De Micheaux, 2010).
#'
#' \bold{Warning:} No missing data is allowed for \code{compResid}, \code{G}, 
#' \code{w} and \code{starResid}.
#'
#' @param compResid a nx1 vector containing the completed residuals
#' @param G a nxs matrix containing the set of SNPs. Each row represents a 
#' different individual and each column represents a separate SNP. The SNP 
#' genotypes should be equal to the number of copies of the minor allele (0, 1 
#' or 2).
#' @param w a sx1 vector of weights for the s SNPs
#' @param ker (default="LIN") Type of kernel matrix: weighted linear ("LIN")
#' or weighted identical-by-state ("IBS")
#' @param asv (default=NULL) Number of approximate eigenvalues to be estimated
#' for the kernel matrix using the implicitly-restarted Lanczos 
#' bidiagonalization implemented in the package \pkg{irlba} (Baglama and 
#' Reichel, 2005). If the spectral decomposition of the matrix is to be 
#' conducted using the R base function \bold{eigen}, \code{asv} can be left as
#' \code{NULL}. This argument has no effect if \code{method} is not equal to "davies".
#' @param method (default="davies") Procedure used to obtain the p-value of the
#' test. "davies" represents the approximation of Davies (1980), "rspMom" 
#' represents the permutation approach based on matching moments described in 
#' Lee et al. (2012), and "rspOrd" represents the standard permutation procedure.
#' 
#' @param starResid (default=NULL) a Bxn \bold{matrix} of permuted residuals
#' used to obtain the p-value of the test following a permutation procedure 
#' (method based on matching moments or standard permutation method). Each row 
#' represents a different permutation sample, and each column represents a 
#' different individual. This argument has no effect if \code{method} is not 
#' equal to "rspOrd" or "rspMom".
#' @param bsw (default=NULL) a vx1 vector containing the lower bounds of the v 
#' sliding windows considered for the SNP-set, taking values between 1 and s
#' @param tsw (default=NULL) a vx1 vector containing the upper bounds of the v
#' sliding windows considered for the SNP-set, taking values between 1 and s
#' @param pos (default=NULL) a sx1 vector of SNP positions
#' @param sf (default=FALSE) logical: indicates whether or not cluster computing
#' is used via the package \pkg{snowfall} in order to reduce wall-clock time. 
#' Initialisation and loading of the package \pkg{gyriq} on all nodes including 
#' master must be called beforehand using the functions \code{sfInit} and 
#' \code{sfLibrary} respectively. See the reference manual of \code{snowfall} 
#' for details. When cluster computing is used, the p-value for each sliding
#' window is computed on a separate node.
#' @param fileOut (default="outGyriq.out") a string containing the name and path
#' of the output file where the results are printed (used only if lower and 
#' upper bounds of sliding windows are also given as input; the file is appended
#' for each sliding window in order to reduce resource wastage)
#' @return If the lower and upper bounds of sliding windows are not provided, 
#' the function produces a list consisting of:
#' @return \item{score}{the score statistic of the test}
#' @return \item{pVal}{the p-value}
#' Otherwise, the function produces a data frame where each row represents a
#' sliding window tested. For each window, the following information is 
#' provided:
#' \itemize{
#' \item \code{FirstSNP}: Rank of the SNP corresponding to the lower bound of
#' the sliding window in the SNP-set
#' \item \code{LastSNP}: Rank of the SNP corresponding to the upper bound of
#' the sliding window in the SNP-set
#' \item \code{winSize}: Number of SNPs in the sliding window
#' \item \code{Start}: Position of the SNP corresponding to the lower bound of
#' the sliding window
#' \item \code{Stop}: Position of the SNP corresponding to the upper bound of
#' the sliding window
#' \item \code{Score}: Score statistic of the association test
#' \item \code{P-value}: P-value of the association test
#' \item \code{Message}: If the calculation of the p-value failed, the 
#' corresponding error message is given. Otherwise, "OK" is displayed.
#' }
#' @author Martin Leclerc <martin.leclerc.5@@ulaval.ca> and Lajmi Lakhal Chaieb 
#' <lakhal@@mat.ulaval.ca>
#' @references Baglama J, Reichel L. 2005. Augmented implicitly restarted 
#' Lanczos bidiagonalization methods. SIAM J Sci Comput 27:19-42.
#' 
#' Davies RB. 1980. The distribution of a linear combination of 
#' \eqn{\chi^2} random variables. J R Stat Soc Ser C 29:323-333.
#' 
#' Lee S, Emond MJ, Bamshad MJ et al. 2012. Optimal unified approach for 
#' rare-variant association testing with application to small-sample 
#' case-control whole-exome sequencing studies. Am J Hum Genet 91:224-237.
#' 
#' Duchesne P, Lafaye De Micheaux P. 2010. Computing the distribution of 
#' quadratic forms: further comparisons between the Liu-Tang-Zhang approximation
#' and exact methods. Comput Stat Data Anal 54:858-862.
#' 
#' Lin X, Zhou Q. 2015. coxKM: Cox Kernel Machine SNP-Set Association Test. R 
#' package version 0.3, URL http://www.hsph.harvard.edu/xlin/software.html#coxkm.
#' 
#' Lin X, Cai T, Wu M, Zhou Q, Liu G, Christiani D and Lin X. 2011. Survival 
#' kernel machine SNP-set analysis for genome-wide association studies. Genetic 
#' Epidemiology 35:620-631.
#' 
#' Cai T, Tonini G and Lin X. 2011. Kernel machine approach to testing the 
#' significance of multiple genetic markers for risk prediction. Biometrics 
#' 67:975-986.
#' @examples
#' data(simGyriq)
#' for (i in seq_along(simGyriq)) assign(names(simGyriq)[i], simGyriq[[i]])
#'
#' cr <- genComplResid(U, Delta, Phi, blkID, m=50, X)
#' testGyriq(cr$compResid, G, w, ker="LIN", asv=NULL, method="davies", 
#' starResid=NULL, bsw, tsw, pos)
#' @importFrom CompQuadForm davies
#' @importFrom irlba irlba
#' @export

testGyriq <- function(compResid, G, w, ker="LIN", asv=NULL, method="davies", 
                      starResid=NULL, bsw=NULL, tsw=NULL, pos=NULL, sf=FALSE, 
                      fileOut="outGyriq.out") {

    if (class(G) != "matrix") stop("G is not a matrix")
    if (sum(is.na(G)) != 0) stop("G cannot have any missing values")
    
    if (sum(apply(G, 2, chkVariation)) == 0) 
        stop("All the SNPs in G have no variation!")
    
    if (sum(is.na(w)) != 0) stop("w cannot have any missing values")
    if (sum(w < 0) != 0) stop("Weights have to be non-negative")
    
    if (length(compResid) != nrow(G)) 
        stop("Dimensions of the completed residuals and G do not match")
    if (length(w) != ncol(G)) 
        stop("Dimensions of w and G do not match")
    
    if (method %in% c("rspMom", "rspOrd")) {
        if (nrow(G) != ncol(starResid))
            stop("Dimensions of G and 'starResid' do not match")        
    }

    if (!(method %in% c("davies", "rspMom", "rspOrd")))
        stop(paste(method, " is not a pre-specified name for 'method'", sep=""))
    
    if (!(ker %in% c("LIN", "IBS")))
        stop(paste(ker, " is not a pre-specified name for 'ker'", sep=""))

    if (ker == "IBS" & method == "davies")
        stop("'IBS' kernel matrix cannot be used with 'davies' method")
                
    if (is.null(bsw) == TRUE) {
        
        W <- diag(w)
        if (ker == "LIN") {
            K <- G %*% W %*% W %*% t(G)
        }
        if (ker == "IBS") {
            given_weight <- 1
            n <- nrow(G)
            p <- ncol(G)
            K <- matrix(rep(0, n * n), nrow = n, ncol = n)
            temp <- .C("Kernel_IBS_Weight", as.integer(as.vector(t(G))), 
                       as.integer(n), as.integer(p), as.integer(given_weight), 
                       as.double(w), as.double(as.vector(K)), 
                       PACKAGE="gyriq")[[6]]
            K <- matrix(temp, nrow=n)
            rm(temp)
        }

        rm(G)
        score <- as.numeric(t(compResid) %*% K %*% compResid)
        resP <- approxPval(score, K, asv, method, starResid)
        return(list(score=signif(score), pVal=signif(resP[[1]])))
        
    } else {
        
        if (length(bsw) != length(tsw)) stop("Dimensions of 'bsw' and 'tsw' do 
                                             not match")
        
        indSW <- as.matrix(seq(along=bsw))
        colNames=c("FirstSNP", "LastSNP", "winSize", "Start", "Stop", "Score",
            "P_value", "Message")
        write(colNames, file=fileOut, ncolumns=11, append=FALSE, sep="\t")
        
        if (sf == TRUE) {
            if ("package:snowfall" %in% search() == FALSE) {
                stop("Package 'snowfall' for cluster computing is not loaded.")
            }
            resList <- snowfall::sfClusterApplyLB(indSW, testSlidWin, compResid,
                                                  G, w, ker, asv, method, 
                                                  starResid, bsw, tsw, pos, 
                                                  fileOut)
            resVec <- unlist(resList)
            nFld <- length(resVec) / length(resList)
            resBM <- matrix(resVec, nrow=length(resList), ncol=nFld, byrow=TRUE)
        } else {
            resBM <- t(apply(indSW, 1, testSlidWin, compResid, G, w, ker, asv, 
                             method, starResid, bsw, tsw, pos, fileOut))
        }
        
        resBM <- as.data.frame(resBM)
        names(resBM) <- c("FirstSNP", "LastSNP", "winSize", "Start", "Stop",
                          "Score", "P_value", "Message")
        return(resBM)
    }
}
