#' a new discriminant Non-Negative Matrix Factorization (dNMF)
#' 
#'  The ndNMF algorithm with the additional Fisher criterion on the cost 
#'  function of conventional NMF was designed to increase class-related
#'  discriminating power.
#'  
#'  This algorithm is based on articles.
#'  \enumerate{
#'  \item Kim, Bo-Kyeong, and Soo-Young Lee. "Spectral Feature Extraction Using dNMF for Emotion Recognition in Vowel Sounds." Neural Information Processing. Springer Berlin Heidelberg, 2013.
#'  \item Lee, Soo-Young, Hyun-Ah Song, and Shun-ichi Amari. "A new discriminant NMF algorithm and its application to the extraction of subtle emotional differences in speech." Cognitive neurodynamics 6.6 (2012): 525-535.
#'  }
#' @param dat a matrix with gene in row and sample in column
#' @param trainlabel the label of sample, like c(1,1,2,2,2)
#' @param r the dimension of expected reduction dimension, with the default value 2
#' @param lambada a relative weighting factor for the discriminant. Default 0.1
#' @param maxIter the maximum iteration of update rules, with the default value 1000
#' @param tol the toleration of coverange, with the default value 1e-7
#' @param log log2 data. Default is TRUE.
#' @param plotit whether plot H (V=WH). Default: FALSE.
#' @param verbose TRUE
#' @param ... to gplots::heatmap.2
#' @author Zhilong Jia and Xiang Zhang
#' @export
#' @import Matrix
#' @examples
#' dat <- rbind(matrix(c(rep(3, 16), rep(8, 24)), ncol=5), 
#' matrix(c(rep(5, 16), rep(5, 24)), ncol=5), 
#' matrix(c(rep(18, 16), rep(7, 24)), ncol=5)) + 
#' matrix(runif(120,-1,1), ncol=5)
#' trainlabel <- c(1,1,2,2,2)
#' 
#' res <- ndNMF(dat, trainlabel, r=2, lambada = 0.1)
#' res$H
#' res$rnk
#' 

ndNMF <- function(dat, trainlabel, r=2, lambada=0.1, maxIter=1000, tol=1e-7, log=TRUE, plotit=FALSE, verbose=FALSE, ...){
    
    eps = .Machine$double.eps
    dat <- as.matrix(dat)
    dat[which(dat==0)] <- 1
    
    if (log){
        dat <- log2(dat + 2) 
    }
    
    nFea = nrow(dat); nSmp = ncol(dat)
    
    # init the H0 and W0 matrix 
    # The 1st row of H is down-regualted genes, 
    # while the 2nd row of H is up-regualted genes)
    H = matrix(runif(r*nSmp, eps), r, nSmp)
    for (i in 1:r){
        H[i,which(trainlabel==names(table(trainlabel))[i])] = H[i,which(trainlabel==names(table(trainlabel))[i])] + sum(H)
    }
    H = pmax(H,0)
    W = matrix(runif(nFea*r, eps), nFea, r)
    # W = W/colSums(W)
    W = pmax(W, eps)
    b = pmax(abs(W %*% H), eps)
    
    obj0 = -10
    lambada = lambada * nFea/r
    
    #averaing matrix Ma over all samples and Mc for wach class
    Ma = matrix(1/nSmp, nrow=nSmp, ncol=nSmp)
    Mc <- list()
    for (i in unique(trainlabel)){
        Mc[[i]] = matrix(1/table(trainlabel)[i], table(trainlabel)[i], table(trainlabel)[i])
    }
    Mc = as.matrix(Matrix::bdiag(Mc))
    
    Hclass = H
    final = Inf
    count = 1
    obj_stack =  vector (mode="numeric", length = maxIter)
    
    while (final > tol && count <= maxIter) {
        
        # update W and H
        W = W * (dat %*% t(H) / (W %*% (H %*% t(H)) + eps) )
        H = H * ( (t(W) %*% dat + lambada * H %*% Mc) / (t(W) %*% W %*% H + lambada * H %*% Ma + eps) )
        if (verbose) print (count)
        
        # H normalization (Formula 11)
        H = as.matrix(H)
        H = H / (rowSums(H)/nSmp + eps)
        
        # Convergence
        b <- pmax(abs( W%*%H ), eps)
        obj1 = norm(dat-b+eps, type="F")^2 - lambada * sum(diag(H %*% (Ma-Mc) %*% t(H)))
        final = abs(obj0-obj1) / (abs(obj0) + eps)
        obj0 = obj1
        if (verbose) print (final)
        obj_stack[count] = obj1
        count = count + 1
        
    }
    
    # to plot H
    if (plotit){
        gplots::heatmap.2(H, scale="row", trace="none", density.info="none", keysize=1, cexCol=0.8, srtCol=30, ...)
    }
    
    list(V=dat, W=W, H=H, rnk=W[,2]-W[,1], trainlabel=trainlabel, count=count,
         final=final, obj_stack=obj_stack, r=r, call=match.call())
}


