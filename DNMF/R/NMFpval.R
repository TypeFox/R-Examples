#' P value for discriminant Non-Negative Matrix Factorization
#' 
#' Estimate the significance of differentially expressed genes in parallel.
#' 
#' P value is caculated based on aatricle, Wang, Hong-Qiang, Chun-Hou Zheng, and Xing-Ming Zhao. 
#' "jNMFMA: a joint non-negative matrix factorization meta-analysis of
#'  transcriptomics data." Bioinformatics (2014): btu679.
#' 
#' @param nmf_res result from DNMF or dNMF
#' @param np number of permutations
#' @param ncores cores used. Default is all the availiable cores
#' @param fdr false discovery rate. Default is FALSE
#' @param top only include top ranked genes. Default is 1000
#' @param verbose verbose
#' @return a matrix with columns rnk, p (and fdr)
#' @export
#' @import parallel
#' @import foreach
#' @import doParallel
#' @author Zhilong Jia
#' @examples
#' dat <- rbind(matrix(c(rep(3, 16), rep(8, 24)), ncol=5), 
#' matrix(c(rep(5, 16), rep(5, 24)), ncol=5), 
#' matrix(c(rep(18, 16), rep(7, 24)), ncol=5)) + 
#' matrix(runif(120,-1,1), ncol=5)
#' trainlabel <- c(1,1,2,2,2)
#' 
#' nmf_res <- ndNMF(dat, trainlabel, r=2, lambada = 0.1)
#' pMat <- NMFpval(nmf_res, np=10, ncores=2, top=4)
#' 
NMFpval <- function (nmf_res, np=100, ncores=parallel::detectCores(), fdr=FALSE, top=1000, verbose=FALSE) {
    
    W <- nmf_res$W
    rownames(W) <- NULL
    rnk0 <- nmf_res$rnk
    rnk <- rnk0[which(abs(rnk0)>sort(abs(rnk0), decreasing=TRUE)[top] )]
    
    
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    if (verbose) {print(paste("getDoParWorkers:", foreach::getDoParWorkers()))}
    strt<-Sys.time()
    
    j <- NULL
    matP <- foreach::foreach(1:np, .combine='cbind') %:%
        foreach::foreach(j=rnk, .combine='c') %dopar% {
            length(which(abs(j)< abs(sample(W[,2], length(rnk0)) - W[,1])))
        }
    print(Sys.time()-strt)
    parallel::stopCluster(cl)
    
    p <- rowMeans(matP)/ length(rnk0)
    if (fdr) {
        fdr <- p.adjust(p, method="fdr")
        dat <- cbind(rnk, p, fdr)
    } else {
        dat <- cbind(rnk, p)
    }
    
    return (dat)
}
