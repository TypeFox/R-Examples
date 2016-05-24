#' Discriminant Non-Negative Matrix Factorization.
#'
#' Discriminant Non-Negative Matrix Factorization, DNMF, is to extend the Non-negative Matrix Factorization algorithm in 
#' order to extract features that enforce not only the spatial locality, but
#'  also the separability between classes in a discriminant manner.
#' 
#' The main algorithm is based on 
#' \href{http://www.ncbi.nlm.nih.gov/pubmed/16722172}{Zafeiriou, S., et al. 
#' (2006) Exploiting discriminant information in 
#' nonnegative matrix factorization with application to frontal face 
#' verification, IEEE transactions on neural networks, 17, 683-695}, 
#' with some \strong{CORRECTIONs}. 
#'
#' @param data a matrix, like expression profilings of some samples. the columns are samples and the rows are gene's expression.
#' @param trainlabel a numeric vector of sample type of all the samples, this vector should ONLY contain 1 and 2 so far and length of it should equal the column (sample) size of data.
#' @param r the dimension of expected reduction dimension, with the default value 2.
#' @param gamma the tradeoff value for the within scatter matrix, with the default value 0.1.  
#' @param delta the tradeoff value for the between scatter matrix, with the default value 1e-4.
#' @param maxIter the maximum iteration of update rules, with the default value 1000.
#' @param log log2 data. Default is TRUE.
#' @param tol the toleration of coverange, with the default value 1e-7.
#' @param plotit whether plot H (V=WH). Default: FALSE.
#' @param checkH whether or not check H. Default: TRUE. This parameter aims to 
#' check whether or not the H safisfy the discriminant metagenes. Usually, this
#' should be TRUE.
#' @param ... to gplots::heatmap.2
#' @import gplots
#' @author Zhilong Jia and Xiang Zhang
#' @export
#' @examples
#' dat <- rbind(matrix(c(rep(3, 16), rep(8, 24)), ncol=5), 
#' matrix(c(rep(5, 16), rep(5, 24)), ncol=5), 
#' matrix(c(rep(18, 16), rep(7, 24)), ncol=5)) + 
#' matrix(runif(120,-1,1), ncol=5)
#' trainlabel <- c(1,1,2,2,2)
#' 
#' DNMF_result <- DNMF(dat, trainlabel, r=2)
#' 
#' 
#' \dontrun{
#' # Gene ranking. dat is the raw read count maatrix with sample in column.
#' 
#' #normalising dat
#' Sizefactors <- DESeq::estimateSizeFactorsForMatrix(dat)
#' dat = sweep(dat, 2, Sizefactors, `/`)
#' 
#' res <- DNMF(dat, trainlabel, r=2)
#' rnk <- res$rnk
#' 
#' #The end of gene ranking exmaples
#' 
#' #Other exmaples
#' DNMF_result <- DNMF(dat, trainlabel, r=2, gamma=0.1, delta=0.0001, plotit=TRUE)
#' }
#' 

DNMF <- function(data,trainlabel, r=2, gamma=0.1, delta=0.0001, maxIter=1000, 
                 tol=1e-7, log=TRUE, plotit=FALSE, checkH=TRUE, ...) {
	
    data <- as.matrix(data)
    data[which(data==0)] <- 1
    
    if (log){
        data <- log2(data + 2) 
    }
    nFea = nrow(data); nSmp = ncol(data)
    eps = .Machine$double.eps
    
    #init the H0 and W0 matrix
    H = matrix(runif(r*nSmp, eps), r, nSmp)
    # The 1st row of H is down-regualted genes, 
    # while the 2nd row of H is up-regualted genes)
    for (i in 1:r){
        H[i,which(trainlabel==names(table(trainlabel))[i])] = H[i,which(trainlabel==names(table(trainlabel))[i])] + sum(H)
    }
    H = pmax(H,0)
    W = matrix(runif(nFea*r, eps), nFea, r)
    W = W/colSums(W)

    #calculate KL divergence of two matrix
    b = pmax(abs(W %*% H), eps)
    obj0 = sum(data*log((data+eps)/(b-eps))-data+b)
    obj1 = obj0 
    
    ##########################
    E = matrix(1, nFea, nSmp)
    # N is just 1/Nr in paper, the weighted matrix
    SmpCount_withinClass = vector(mode="numeric", length(trainlabel))
    SmpCount_withinClass = as.vector(table(trainlabel)[trainlabel])
    N = matrix(1/SmpCount_withinClass, r, nSmp, byrow=T)     

    final = Inf
    count = 1
    Hclass = H  
    obj_stack =  vector (mode="numeric", length = maxIter)
    
while (final > tol && count <= maxIter) {

    #update H with the objective function includes KL divergence

    for(i in unique(trainlabel)){
        Hclass[,which(trainlabel==i)] = matrix( rep(rowSums(H[,which(trainlabel==i)]), length(which(trainlabel==i))), r, length(which(trainlabel==i)))
    }
    Hclass = Hclass - H
    Hsum = matrix(rep(rowSums(H),ncol(H)),r, ncol(H)) - H

    tmp_a = 4*gamma + 4*(delta/nSmp - (gamma+delta)*N)   #2a
    tmp_b = 1 + 2*delta/nSmp * Hsum - 2*(delta+gamma)*N * Hclass
    tmp_c = - t(W) %*% (data / ( W %*% H + eps) ) * H
    H = (sqrt(tmp_b^2 - 2*tmp_a*tmp_c) - tmp_b ) /tmp_a
    H = pmax(H, eps)  

    #######################################
    #update W
    W = (W/(E%*%t(H)))*(data/(W%*%H)) %*% t(H)
    H = diag(colSums(W)) %*% H
    W = W/ matrix(rep(colSums(W), each=nrow(W)), nrow(W), ncol(W))
    W = pmax(W,eps)
    
    obj2 = obj1
    b = pmax(abs(W%*%H),eps)
    obj1 = sum( data*log((data+eps)/(b-eps)) - data + b )
    final = abs(obj1-obj2) / abs(obj1-obj0)
    #obj_stack[count] = final
    obj_stack[count] = obj1
    count = count + 1
}
    # to plot H
    if (plotit){
        gplots::heatmap.2(H, scale="row", trace="none", density.info="none", keysize=1, cexCol=0.8, srtCol=30, ...)
    }
    

    # check H. Setting down-regulted metagene in row 1 of H, 
    # while up-regulated metagene in row 2 of H.
    l1 <- apply(H[,which(trainlabel == names(table(trainlabel))[1])], 1, mean)
    l2 <- apply(H[,which(trainlabel == names(table(trainlabel))[2])], 1, mean)
    meanH <- cbind(l1, l2)
    if (checkH){
        if (l1[1] < l2[1] && l1[2] > l2[2]) {
            H <- rbind(H[2,], H[1,])
            W1 <- cbind(W[,2], W[,1]) 
        } else if ((l1[1] < l2[1] && l1[2] < l2[2]) || (l1[1] > l2[1] && l1[2] > l2[2])) {
            stop("Failed. Run DNMF again after restart R.")
        }
    }
    
    list(V=data, W=W, H=H, rnk=W[,2]-W[,1], trainlabel=trainlabel, delta=delta, gamma=gamma, count=count,
         final=final, obj_stack=obj_stack, r=r, call=match.call())
}


