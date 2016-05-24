#'Eigenvalue plot
#' 
#' @details
#' Plots eigenvalues of the correlation matrix and overlays the Marchenko-Pastur
#' density on top of it. There is a shap cutoff for the density. We are concerned
#' with eigenvalues beyond this cutoff. Paramters used for plotting are added 
#' to the plot
#' 
#' @param x model of the type RMT obtained by fitting an RMT model to the data
#' @param y unused
#' @param ... additional arguments unused
#' @author Rohit Arora
#' @examples 
#' \dontrun{
#'  data("largereturn")
#'  model <- estRMT(largesymdata)
#'  plot(model)
#' }
#' 
#' @method plot RMT
#' @export
#' 
plot.RMT <- function(x, y, ...){
    
    lambdas <- x$eigVals; Q <- x$Q; sigma.sq <- x$var 
    lambda.max <- x$lambdascutoff 
    
    p <- ggplot(data=data.frame(lambdas)) + 
        geom_histogram( aes_string(x = 'lambdas', y='..density..'),
                        breaks=seq(min(lambdas)-1,1+max(lambdas),0.5), 
                        colour="black", fill="white") +
        stat_function(fun = dmp, args=list(svr = Q, var=sigma.sq), 
                      aes(colour = 'MP density')) + xlab("Eigenvalues") +
        labs(title="Actual vs Fitted Marchenko-Pastur") + ylim(0,1.5) + 
        theme(plot.title = element_text(size = 20, face = "bold", vjust = 1),
              axis.title=element_text(size=14,face="bold")) + 
        annotate('text', x = 10, y = 0.9, 
                 label = paste("sigma^{2} == ", round(sigma.sq,3)), parse=TRUE) +
        annotate('text', x = 10, y = 1, 
                 label = paste("Q == ", round(Q,3)), parse=TRUE) + 
        annotate('text', x = 10, y = 0.78, 
                 label = paste("lambda[max] ==", round(lambda.max,3)), parse=TRUE) + 
        scale_colour_manual("", values = c("red"))
    
    options(warn = -1)
    print(p)
    options(warn = 0)
    p
}

#' Denoising of Covariance matrix using Random Matrix Theory
#' 
#' @details
#' This method takes in data as a matrix or an xts object. It then
#' fits a marchenko pastur density to eigenvalues of the correlation matrix. All
#' eigenvalues above the cutoff are retained and ones below the cutoff are
#' replaced such that the trace of the correlation matrix is 1 or non-significant
#' eigenvalues are deleted and diagonal of correlation matrix is changed to 1. 
#' Finally, correlation matrix is converted to covariance matrix.
#' 
#' @importFrom Matrix nearPD
#' @importFrom RMTstat dmp qmp
#' @importFrom foreach "%dopar%" foreach registerDoSEQ
#' @importFrom parallel detectCores stopCluster makeCluster clusterEvalQ clusterExport
#' @importFrom doParallel registerDoParallel
#' 
#' @param  R xts or matrix of asset returns
#' @param  Q ratio of rows/size. Can be supplied externally or fit using data
#' @param  cutoff takes two values max/each. If cutoff is max, Q is fitted and 
#'          cutoff for eigenvalues is calculated. If cutoff is each, Q is set to
#'          row/size. Individual cutoff for each eigenvalue is calculated and used
#'          for filteration. 
#' @param eigenTreat takes 2 values, average/delete. If average then the noisy 
#'        eigenvalues are averged and each value is replaced by average. If delete
#'        then noisy eigenvalues are ignored and the diagonal entries of the 
#'        correlation matrix are replaced with 1 to make the matrix psd.
#' @param numEig number of eigenvalues that are known for variance calculation.
#'        Default is set to 1. If numEig = 0 then variance is assumed to be 1.
#' @param parallel boolean to use all cores of a machine.
#' @examples 
#' \dontrun{
#'  data("largereturn")
#'  model <- estRMT(largesymdata, numEig = 0)  
#' }        
#'        
#' @author Rohit Arora
#' 
#' @export
#' 
estRMT <- function(R, Q =NA, cutoff = c("max", "each"), 
                    eigenTreat = c("average", "delete") , numEig=1, 
                   parallel = TRUE) {
    .data <- if(is.xts(R)) coredata(R) else as.matrix(R)
    T <- nrow(.data); M <- ncol(.data) 
    if (T < M) stop("Does not work when T < M")
    
    if(!is.na(Q)) if(Q < 1) stop("Does not work for Q<1")
    
    cutoff <- cutoff[1]; if(!cutoff %in% c("max", "each")) stop("Invalid cutoff")
    if(cutoff == "each") Q <- T/M
    
    eigenTreat <- eigenTreat[1]; 
    if(!eigenTreat %in% c("average", "delete")) stop("Invalid eigenTreat option")
    
    if (numEig < 0) stop("Number of eigenvalues must be non-negative")
    
    #eigenvalues can be negative. To avoid this e need a positive-definite matrix 
    S <- cov(.data); S <- as.matrix(nearPD(S)$mat)
    D <- diag(diag(S)); C <- cov2cor(S); 
    
    # Marchenko Pastur density is defined for eigenvalues of correlation matrix
    eigen.C <- eigen(C,symmetric=T)
    lambdas <- eigen.C$values; sigma.sq <- mean(lambdas)
    
    #minimize log-likelihood. 
    loglik.marpas <- function(theta, sigma.sq) {
        
        Q <- theta
        val <- sapply(lambdas,     
                      function(x) dmp(x,svr = Q, var=sigma.sq))
        
        val <- val[val > 0]
        ifelse(is.infinite(-sum(log(val))), .Machine$double.xmax, -sum(log(val)))        
    }
    
    sigma.sq <- 1 - sum(head(lambdas,numEig))/M
    
    if( is.na(Q) && cutoff != "each") {
        lb <- 1; ub <- max(T/M,5)
        if(parallel) {
          cl <- makeCluster(detectCores())
          registerDoParallel(cl)
          clusterEvalQ(cl, library(RMTstat))
        }

        '%exectype%' <- if (parallel) get('%dopar%') else get('%do%')
        
        starts <- seq(lb, ub, length.out = 50)
        fit.marpas <- foreach(start = starts, .combine = rbind) %exectype% 
            optim(par = start, fn = loglik.marpas, method = "L-BFGS-B", 
                  lower = lb, upper = ub, sigma.sq = sigma.sq)   
        
        if(parallel) stopCluster(cl)
        
        idx <- grep("CONVERGENCE",unlist(fit.marpas[,"message"]))
        vals <- fit.marpas[idx,c("par","value")]
        Q <- unlist(vals[which.min(vals[,"value"]),"par"])    
    }

    lambda.max <- qmp(1, svr=Q, var = sigma.sq)  
    # now that we have a fit. lets denoise eigenvalues below the cutoff
    
    idx <- if(cutoff == "max") 
        which(lambdas > lambda.max)
    else if(cutoff == "each")
    {
        cutoff.each <- sapply(2:length(lambdas), function(i) {
            eigr <- lambdas[i:M]
            mean(eigr)*(1 + (M - i + 1)/T + 2*sqrt((M - i + 1)/T))
        })
        
        c(1, 1 + which(lambdas[-1] > cutoff.each))    
    }

    if (length(idx) == 0) return(S)
    
    val <- eigen.C$values[idx]; vec <- eigen.C$vectors[,idx,drop=FALSE]
    sum <- 0; for (i in 1:ncol(vec)) sum <- sum + val[i]*vec[,i] %*% t(vec[,i])
    
    # trace of correlation matrix is 1. Use this to determine all the remaining
    # eigenvalues
    
    lambdas.cleaned <- c()
    clean.C <- if (eigenTreat == "average") {
        lambdas.cleaned <- c(val, rep(1,M))
        sum + sum(eigen.C$values[-idx])/M * diag(rep(1,M))
    } else if (eigenTreat == "delete") {
        lambdas.cleaned <- c(val, rep(0,M))
        diag(sum) <- 1
        sum
    }
    
    # convert correlation to covariance matrix and return
    clean.S <- D^0.5 %*% clean.C %*% D^0.5
    fit <- list(cov = clean.S, Q = Q, var = sigma.sq, eigVals = lambdas, 
                eigVals.cleaned = lambdas.cleaned, lambdascutoff = lambda.max)
    
    class(fit) <- "RMT"
    fit
}
