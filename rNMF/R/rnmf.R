#' Robust Penalized Nonnegative Matrix Factorization (rNMF).
#' 
#' rNMF performs robust penalized nonnegative matrix factorizaion on a nonnegative matrix X to obtain low dimensional nonnegative matrices W and H, such that X ~ WH, while detecting and trimming different types of outliers in X.
#'
#' rNMF decomposes a nonnegative p by n data matrix X as X ~ WH and detect and trims outliers. Here W and H are p by k and k by n nonnegative matrices, respectively; and k <= min\{p, n\} is the dimension of the subspace to which X is projected. The objective function is
#' 
#' ||X - WH||_{gamma} + alpha * ||W||_2^2 + beta * sum(|H_.j|)^2
#' 
#' where alpha controls the magnitude of W, and beta controls the sparsity of H. The algorithm iteratively updates W, H and the outlier set zeta with alternating conditional nonnegative least square fittings until convergence.
#'
#' Four variations of trimming are included in the algorithm: "cell", "row", "column" and "smooth". Specifically, the "cell" variation trims individual cell-wise outliers while "row" and "column" variations trim entire row or column outliers. The fourth variation "smooth" fills the cells that are declared outliers in each iteration by the average of the surrounding cells. 
#' 
#' @param X a 'p by n' nonnegative numerical matrix to be decomposed into X ~ WH.
#' @param k the integer dimension to which X is projected (equals the number of columns of W and the number of rows of H, and smaller than min(p,n)). Default = 5.  
#' @param alpha a nonnegative number that controls the magnitude of W. The larger the alpha is, the higher penalty on ||W|| is, or ||W|| is forced to be smaller. Default = 0, representing no penalty on the magnitude ||W||.
#' @param beta a nonnegative number that controls the sparsity of H. Default= 0. The larger the beta is, the higher penalty on sum_i ||H_j|| is, where H_j is the j-th column of H. Default = 0, representing no penalty on the sparsity of H. 
#' @param maxit the maximum number of iterations. Default = 20. The algorithm is done as follows: 1) fits W given H, then trims outliers in X, i.e. ones with large residuals from X - WH, then refits W without the outliers [this step can be repeated  'nreg' times, currently nreg = 1],  then 2) repeat as in 1) with the roles of W and H swapped,  and then iterate through 1) and 2) until convergence. Default = 20, allowing a maximum of 20 pairs of 1) and 2).
#' @param tol the convergence tolerance. We suggest it to be a positive number smaller than 0.01. Default = 0.001.
#' @param gamma a trimming percentage in (0, 1) of X. Default = 0.05 will trim 5\% of elements of X (measured by cell,  row or column as specified in 'variation'). If trim=0,  there is no trim; so rNMF then performs the regular NMF.
#' @param ini.W the initialization of the left matrix W, which is a "p by k" nonnegative numeric matrix. Default = NULL directs the algorithm to randomly generate an ini.W.
#' @param ini.zeta a "p by n" logical matrix of True or False, indicating the initial locations of the outliers. The number of "TRUE"s in ini.zeta  must be less than or equal to m = the rounded integer of (gamma * p * n).  Default = NULL, initializes the cells as TRUE with the m largest residuals after the first iteration. Required only for "cell" trimming option.  
#' @param my.seed the random seed for initialization of W. If left to be NULL(default) a random seed is used.  Required only if ini.W is not NULL.
#' @param variation a character string indicating which trimming variation is used. The options are: 'cell' (default), 'col', 'row' and 'smooth'.
#' @param quiet default = FALSE indicating a report would be given at the end of execution. If quiet = TRUE, no report is provide at the end.
#' @param nreg the number of inner loop iterations [see 'maxit' above] to  find outliers in X, given either H or W. Default = 1, is currently only implemented option in the "cell" variation.
#' @param showprogress default = TRUE, shows a progress bar during iterations. If showprogress = FALSE, no progress bar is shown.
#' 
#' @return An object of class 'rnmf', which is a list of the following items:
#' \itemize{
#' \item W: left matrix of the decomposition X ~ WH, columns of which (i.e. W) are basis vectors of the low dimension projection.
#' \item H: right matrix of the decomposition X ~ WH, columns of which (i.e. W) are low dimensional encoding of the data.
#' \item fit: the fitted matrix W \%*\% H.
#' \item trimmed: a list of locations of trimmed cells in each iteration.
#' \item niter: the number of iterations performed.
#' }
#'
#' @export
#' 
#' @examples
#' ## Load a clean single simulated tumor image.
#' data("tumor")
#' image(tumor) # look at the tumor image
#' dim(tumor) # it is a '70 by 70' matrix
#' ## Add 5\% corruptions.
#' tumor.corrupted <- tumor
#' set.seed(1)
#' tumor.corrupted[sample(1:4900, round(0.05 * 4900), replace = FALSE)] <- 1
#' ## Run rnmf with different settings
#' # No trimming
#' res.rnmf1 <- rnmf(X = tumor.corrupted, gamma = FALSE, my.seed = 1)
#' # 6 percent trimming, low dimension k = 5 (default)
#' res.rnmf2 <- rnmf(X = tumor.corrupted, tol = 0.001, gamma = 0.06, my.seed = 1)
#' # add sparsity constraint of H (beta = 0.1) with k = 10, and the "smooth" variation.
#' res.rnmf3 <- rnmf(X = tumor.corrupted, k = 10, beta = 0.1,
#'                   tol = 0.001, gamma = 0.06, my.seed = 1, 
#'                   variation = "smooth", maxit = 30)
#'
#' ## Show results:
#' par(mfrow = c(2,2), mar = c(2,2,2,2))
#' image(tumor.corrupted, main = "Corrupted tumor image", xaxt = "n", yaxt = "n")
#' image(res.rnmf1$fit, main = "rnmf (no trimming) fit", xaxt = "n", yaxt = "n")
#' image(res.rnmf2$fit, main = "rnmf (cell) fit 2", xaxt = "n", yaxt = "n")
#' image(res.rnmf3$fit, main = "rnmf (smooth) fit 3", xaxt = "n", yaxt = "n")

rnmf <- function(X, k = 5, alpha = 0, beta = 0, maxit = 20, tol = 0.001,
    gamma = FALSE, ini.W = NULL, ini.zeta = NULL, my.seed = NULL, 
    variation = "cell", quiet = FALSE, nreg = 1, showprogress = TRUE)
{
    tic <- proc.time() # Start a clock.
    p <- nrow(X)
    n <- ncol(X)
    ## Checks if the arguments are valid. Display an error message if not.
    checkargs(X, k, alpha, beta, maxit, tol, gamma, ini.W, ini.zeta, my.seed, 
              variation, quiet, nreg, p, n)
    
    ## Create a data frame of 4 columns: value = entries of X; 
    ## (x,y) = coordinates; outs = is it an outlier?
    X.f <- data.frame(value = as.vector(X), x = rep(1 : p, n),
                      y = rep(1 : n, each = p), outs = FALSE) 

    ## to.trim (list) stores entries trimmed in each iteration.
    if(gamma > 0) to.trim <- vector("list")
    else to.trim <- NULL
    ## 'obj' (vector) stores values of the objective function in each iteration.
    obj <- 1

    ## Initialize W
    if(missing(ini.W)){
        if(missing(my.seed)){
            W <- initM(large = max(X), nrow = p, ncol = k, small = 0)
        }else{ 
            W <- initM(large = max(X), nrow = p, ncol = k, small = 0,
                       my.seed = my.seed)
        }
    }else{
        W <- ini.W
    }
    ## Creates a progress bar.
    if(showprogress) pb <- txtProgressBar(min = 0, max = maxit, style = 3)
    
    if(variation == "cell" & gamma > 0){
      ##--------------------------
      ## Cell trimming
      ##--------------------------
        ## Initialize zeta (logic matrix the same size as X)
        if(!is.null(ini.zeta)){
            if(sum(ini.zeta) < round(gamma * p * n)) {
                warning("The percentage of FALSES (outliers) in ini.zeta 
                        is smaller than 'gamma'; Other outliers are randomly picked.")
                more.outliers <- sample(which(!ini.zeta), 
                                    round(gamma * p * n) - sum(ini.zeta))
                ini.zeta[more.outliers] <- TRUE
            }else{}
            zeta <- ini.zeta
        }else{
            zeta <- matrix(FALSE, nrow = p, ncol = n)
            zeta[sample(1 : (p * n), round(gamma * p * n))] <- TRUE ## randomize zeta
        }

        ## Start iterations.
        for(i in 1 : maxit){
            if(showprogress) setTxtProgressBar(pb, i) ## update the progress bar
            ##------Stage 1------##
            ## Given W, fit H
            H <- Nnls.trimH(W, X, !zeta, beta, k, n)
            for(j in 1:nreg){
                ## Find residuals
                R <- abs(X - W %*% H)
                to.trim[[i]] <- order(R, decreasing = TRUE)[1 : round(gamma * n * p)]
                ## to.trim[[i]] = which(rank(R) > round((1 - gamma) * n * p))
                ## Update zeta
                zeta <- matrix(FALSE, nrow = p, ncol = n)
                zeta[to.trim[[i]]] <- TRUE
                ## Refit H
                H <- Nnls.trimH(W, X, !zeta, beta, k, n)
            }
            
            ##------Stage 2------##
            ## Given H, Fit W
            W <- Nnls.trimW(H, X, !zeta, alpha, k, p)
            for(j in 1:nreg) {
                ## Find residuals
                R <- abs(X - W %*% H)
                ##to.trim[[i]] = which(rank(R) > round((1 - gamma) * n * p))
                to.trim[[i]] <- order(R, decreasing=TRUE)[1 : round(gamma * n * p)]
                ## Update zeta
                zeta <- matrix(FALSE, nrow = p, ncol = n)
                zeta[to.trim[[i]]] <- TRUE
                ## Refit W
                W <- Nnls.trimW(H, X, !zeta, alpha, k, p)
                J <- nmlz(W) # Find the normalizing matrix of W
                W <- W %*% J; H = solve(J) %*% H
            }
            obj[i] <- l2((X - W %*% H)[!zeta]) + l2(W) + sum(colSums(abs(H))^2)

            ## Check convergence
#             if(i > 1){
#                 if(all(to.trim[[i]] == to.trim[[i - 1]]) &
#                    sum((W - W.prev)^2) / sum(W.prev^2) < tol) break
#             }else{}
            if(i > 1){
              if(sum((W - W.prev)^2) / sum(W.prev^2) < tol) break
            }else{}
            W.prev <- W
        }
    }else{
      ##--------------------------
      ## Row, Col and Smooth trimming
      ##--------------------------  
        for(i in 1 : maxit){
            flag <- FALSE
            if(showprogress) setTxtProgressBar(pb, i)  # update the progress bar
            ##------Stage 1------##
            #Given W, find H in ||MH - C||
            C <- rbind(X, matrix(0,1,n))
            M <- rbind(W, sqrt(beta) * matrix(1, 1, k))
            H <- Nnls(M,C)
            if(gamma > 0){ ## Trim
                R <- (X - W %*% H)^2
                if(variation == "col"){
                    # No change here
                }else if(variation == "smooth"){
                    to.trim[[i]] <- which(rank(abs(R)) > round((1 - gamma) * n * p))
                    X.s <- matrix(smoothing2(X.f, to.trim[[i]], p,n, frame = TRUE)$value,p,n)
                    ## Fit H with trimmed X and W. 
                    C.s <- rbind(X.s, matrix(0, nrow = 1, ncol = n))
                    H.prev <- H
                    H <- Nnls(M,C.s)
                }else if(variation == "row"){
                    rsqr <- apply(R, 1, l2) ## Find SS of residuals of rows
                    to.trim[[i]] <- which(rank(rsqr) > round((1 - gamma) * p))
                    X.trim <- X[-to.trim[[i]],]
                    W.trim <- W[-to.trim[[i]],]
                    ## Fit H with trimmed X and W. 
                    M.trim <- rbind(W.trim, sqrt(beta) * matrix(1, nrow = 1, ncol = k))
                    C.trim <- rbind(X.trim, matrix(0, nrow = 1, ncol = n))
                    H.prev <- H
                    H <- Nnls(M.trim, C.trim)
                }else{
                    stop("Wrong mode. Try one of the following: 'cell', 'col', 'row', 'smooth'")
                }
            }else{}
            
            ##------Stage 2------##
            ##Given H, find W in min||MW^T - C||.
            C <- rbind(t(X), matrix(0, nrow = k, ncol = p))
            M <- rbind(t(H), sqrt(alpha) * diag(k))
            W <- t(Nnls(M, C))
            if(gamma > 0){
                R <- (X - W %*% H)^2
                if(variation == "col"){
                    rsqr <- apply(R, 2, l2) ## Find SS of residuals of columns
                    to.trim[[i]] <- which(rank(rsqr) > round((1 - gamma) * n))
                    X.trim <- X[,-to.trim[[i]]]
                    H.trim <- H[,-to.trim[[i]]]
                    ## Fit W with trimmed X and H. 
                    M.trim <- rbind(t(H.trim), sqrt(alpha) * diag(k))
                    C.trim <- rbind(t(X.trim), matrix(0, nrow = k, ncol = p))
                    W <- t(Nnls(M.trim,C.trim))
                }else if(variation == "smooth"){
                    R <- (X - W %*% H)^2
                    to.trim[[i]] <- which(rank(abs(R)) > round((1 - gamma) * n * p))
                    X.s <- matrix(smoothing2(X.f, to.trim[[i]], p, n, frame = TRUE)$value,p,n)
                    ## Fit W with trimmed X and H. 
                    C.s <- rbind(t(X.s), matrix(0, nrow = k, ncol = p))
                    W <- t(Nnls(M,C.s))
                }else if(variation == "row"){
                    ## No change
                }else{
                    stop("Wrong mode. Try one of the following: 'cell', 'col', 'row', 'smooth'")
                }
            }
            ##J = nmlz(W) # Find the normalizing matrix of W
            ##W = W %*% J; H = solve(J) %*% H
            
            ## Convergence?
            if(i > 1){
                if(setequal(to.trim[[i]], to.trim[[i-1]])){
                    if(sum((W - W.prev) ^ 2) / sum(W.prev ^ 2) < tol) break
                }
            }else{}
            W.prev <- W
        }
    }
    ## Close the link to the progress bar. 
    if(showprogress) close(pb)
    fit <- W %*% H
    
    ## Print report
    if(!quiet){
        if(!gamma){
            cat("Done. Time used:","\n")
            print(proc.time() - tic)
            cat("No trimming.\n",
                "Input matrix dimension: ",p," by ",n, "\n",
                "Left matrix: ", p," by ", k,". Right matrix: ", k," by ", n,"\n",
                "alpha = ", alpha,". beta = ",beta,"\n",
                "Number of max iterations = ",maxit,"\n",
                "Number of iterations = ", i, "\n", sep = ""
                )
        }else{
            cat("Done. Time used: ","\n")
            print(proc.time() - tic)
            cat("\n Trimming mode = \"", variation, "\". Proportion trimmed: ",gamma, "\n",
                "Input matrix dimension: ",p," by ",n, "\n",
                "Left matrix: ",p," by ",k,". Right matrix: ",k," by ",n,"\n",
                "alpha = ",alpha,". beta = ",beta,"\n",
                "Number of max iterations = ",maxit,"\n",
                "Number of iterations = ",i,"\n", sep = ""
                )
        }
    }
    return(invisible(list(W = W, H = H, fit = fit,
                          trimmed = to.trim, niter = i)))
}

