#' SGCCA extends RGCCA to address the issue of variable selection. Specifically, 
#' RGCCA is combined with an L1-penalty that gives rise to Sparse GCCA (SGCCA) which 
#' is implemented in the function sgcca().
#' Given \eqn{J} matrices \eqn{\mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}}, that represent 
#' \eqn{J} sets of variables observed on the same set of \eqn{n} individuals. The matrices 
#' \eqn{\mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}} must have the same number of rows, but may 
#' (and usually will) have different numbers of columns. Bblocks are not necessarily fully connected
#' within the SGCCA framework. Hence the use of SGCCA requires 
#' the construction (user specified) of a design matrix (\eqn{\mathbf{C}}) that characterizes 
#' the connections between blocks. Elements of the (symmetric) design matrix \eqn{\mathbf{C} = (c_{jk})} 
#' areequal to 1 if block \eqn{j} and block \eqn{k} are connected, and 0 otherwise. 
#' Hence, the use of SGCCA requires the construction 
#' (user specified) of a design matrix (\eqn{\mathbf{C}}) which characterizes 
#' the connections between blocks. The SGCCA algorithm is very similar to the RGCCA algorithm and keeps the same monotone 
#' convergence properties (i.e. the bounded criteria to be maximized increases 
#' at each step of the iterative procedure).
#' Moreover, using a deflation strategy, sgcca() enables computation of several SGCCA block 
#' components (specified by ncomp) for each block. Block components for each block are guaranteed to be orthogonal 
#' when using this deflation strategy. The so-called symmetric deflation is considered in this implementation,
#' i.e. each block is deflated with respect to its own component. 
#' Moreover, we stress that the numbers of components per block could differ 
#' from one block to another. 
#' @param A  A list that contains the \eqn{J} blocks of variables \eqn{\mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}}.
#' @param C  A design matrix that describes the relationships between blocks (default: complete design).
#' @param c1 Either a \eqn{1 \times J} vector or a \eqn{\max (ncomp) \times J} matrix encoding the L1 constraints applied to the outer weight vectors. Elements of c1 vary between 0 and 1 (larger values of c1 correspond to less penalization). 
#' If c1 is a vector, L1-penalties are the same for all the weights corresponding to the same block but different components: 
#' \deqn{\forall h, \|a_{j,h}\|_{\ell_1} \leq c_1[j] \sqrt{p_j},}
#' with \eqn{p_j} the number of variables of \eqn{\mathbf{X}_j}.
#' If c1 is a matrix, each row \eqn{h} defines the constraints applied to the weights corresponding to components \eqn{h}:
#' \deqn{\forall h, \|a_{j,h}\|_{\ell_1} \leq c_1[h,j] \sqrt{p_j}.}
#' @param ncomp  A \eqn{1 \times J} vector that contains the numbers of components for each block (default: rep(1, length(A)), which means one component per block).
#' @param scheme Either  "horst", "factorial" or "centroid" (Default: "centroid").
#' @param scale	If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
#' @param init Mode of initialization use in the SGCCA algorithm, either by Singular Value Decompostion ("svd") or random ("random") (default : "svd").
#' @param bias A logical value for biaised or unbiaised estimator of the var/cov.
#' @param verbose  Will report progress while computing if verbose = TRUE (default: TRUE).
#' @param tol Stopping value for convergence.
#' @return \item{Y}{A list of \eqn{J} elements. Each element of Y is a matrix that contains the SGCCA components for each block.}
#' @return \item{a}{A list of \eqn{J} elements. Each element of a is a matrix that contains the outer weight vectors for each block.}
#' @return \item{astar}{A list of \eqn{J} elements. Each element of astar is a matrix defined as Y[[j]][, h] = A[[j]]\%*\%astar[[j]][, h]}
#' @return \item{C}{A design matrix that describes the relationships between blocks (user specified).}
#' @return \item{scheme}{The scheme chosen by the user (user specified).}
#' @return \item{c1}{A vector or matrix that contains the value of c1 applied to each block \eqn{\mathbf{X}_j}, \eqn{ j=1, \ldots, J} and each dimension (user specified).}
#' @return \item{ncomp}{A \eqn{1 \times J} vector that contains the number of components for each block (user specified).}
#' @return \item{crit}{A vector that contains the values of the objective function at each iterations.}
#' @return \item{AVE}{Indicators of model quality based on the Average Variance Explained (AVE): AVE(for one block), AVE(outer model), AVE(inner model).}
#' @references Tenenhaus et al. Variable Selection For Generalized Canonical Correlation Analysis. 2013. Submitted to Biostatistics.
#' @title Variable Selection For Generalized Canonical Correlation Analysis (SGCCA)
#' @examples
#'
#' #############
#' # Example 1 #
#' #############
#' \dontrun{
#' # Download the dataset's package at http://biodev.cea.fr/sgcca/.
#' # --> gliomaData_0.4.tar.gz
#' 
#' require(gliomaData)
#' data(ge_cgh_locIGR)
#' 
#' A <- ge_cgh_locIGR$multiblocks
#' Loc <- factor(ge_cgh_locIGR$y) ; levels(Loc) <- colnames(ge_cgh_locIGR$multiblocks$y)
#' C <-  matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
#' tau = c(1, 1, 0)
#' 
#' # rgcca algorithm using the dual formulation for X1 and X2 
#' # and the dual formulation for X3
#' A[[3]] = A[[3]][, -3]
#' result.rgcca = rgcca(A, C, tau, ncomp = c(2, 2, 1), scheme = "factorial", verbose = FALSE)
#' # sgcca algorithm
#' result.sgcca = sgcca(A, C, c1 = c(.071,.2, 1), ncomp = c(2, 2, 1), 
#'                      scheme = "centroid", verbose = FALSE)
#' 
#' ############################
#' # plot(y1, y2) for (RGCCA) #
#' ############################
#' layout(t(1:2))
#' plot(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[2]][, 1], col = "white", xlab = "Y1 (GE)", 
#'      ylab = "Y2 (CGH)", main = "Factorial plan of RGCCA")
#' text(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[2]][, 1], Loc, col = as.numeric(Loc), cex = .6)
#' plot(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[1]][, 2], col = "white", xlab = "Y1 (GE)", 
#'      ylab = "Y2 (GE)", main = "Factorial plan of RGCCA")
#' text(result.rgcca$Y[[1]][, 1], result.rgcca$Y[[1]][, 2], Loc, col = as.numeric(Loc), cex = .6)
#' 
#' ############################
#' # plot(y1, y2) for (SGCCA) #
#' ############################
#' layout(t(1:2))
#' plot(result.sgcca$Y[[1]][, 1], result.sgcca$Y[[2]][, 1], col = "white", xlab = "Y1 (GE)", 
#'      ylab = "Y2 (CGH)", main = "Factorial plan of SGCCA")
#' text(result.sgcca$Y[[1]][, 1], result.sgcca$Y[[2]][, 1], Loc, col = as.numeric(Loc), cex = .6)
#' 
#' plot(result.sgcca$Y[[1]][, 1], result.sgcca$Y[[1]][, 2], col = "white", xlab = "Y1 (GE)", 
#'      ylab = "Y2 (GE)", main = "Factorial plan of SGCCA")
#' text(result.sgcca$Y[[1]][, 1], result.sgcca$Y[[1]][, 2], Loc, col = as.numeric(Loc), cex = .6)
#' 
#' # sgcca algorithm with multiple components and different L1 penalties for each components 
#' # (-> c1 is a matrix)
#' init = "random"
#' result.sgcca = sgcca(A, C, c1 = matrix(c(.071,.2, 1, 0.06, 0.15, 1), nrow = 2, byrow = TRUE), 
#'                      ncomp = c(2, 2, 1), scheme = "factorial", scale = TRUE, bias = TRUE, 
#'                      init = init, verbose = FALSE)
#' # number of non zero elements per dimension
#' apply(result.sgcca$a[[1]], 2, function(x) sum(x!=0)) 
#'      #(-> 145 non zero elements for a11 and 107 non zero elements for a12)
#' apply(result.sgcca$a[[2]], 2, function(x) sum(x!=0)) 
#'      #(-> 85 non zero elements for a21 and 52 non zero elements for a22)
#' init = "svd"
#' result.sgcca = sgcca(A, C, c1 = matrix(c(.071,.2, 1, 0.06, 0.15, 1), nrow = 2, byrow = TRUE), 
#'                      ncomp = c(2, 2, 1), scheme = "factorial", scale = TRUE, bias = TRUE, 
#'                      init = init, verbose = FALSE)}
#'@export sgcca


sgcca <- function (A, C = 1-diag(length(A)), c1 = rep(1, length(A)), ncomp = rep(1, length(A)), scheme = "centroid", scale = TRUE, init = "svd", bias = TRUE, tol = .Machine$double.eps, verbose = FALSE){
    if ( any(ncomp < 1) ) stop("One must compute at least one component per block!")
    if ( any((c1<0) | (c1>1)) ) stop("c1 must vary between 0 and 1!")
    
    #-------------------------------------------------------
    if ((scheme != "horst" ) & (scheme != "factorial") & (scheme != "centroid")) {
      stop("ERROR : choose one of the three following schemes: horst, centroid or factorial")
    }else{
      if (verbose) cat("Computation of the SGCCA block components based on the",scheme,"scheme \n")
    }
    #-------------------------------------------------------
    
    if (scale == TRUE) A = lapply(A, function(x) scale2(x, bias = bias))
    
    ####################################
    # sgcca with 1 component per block #
    ####################################
          
    # ndefl number of deflation per block
    ndefl <- ncomp-1
    N <- max(ndefl)
    J <- length(A)
    pjs <- sapply(A,NCOL)
    nb_ind <- NROW(A[[1]])
    AVE_X = list()
    AVE_outer <- rep(NA,max(ncomp))
    if (N == 0) {
        result <- sgccak(A, C, c1, scheme, init = init, bias = bias, tol = tol, verbose = verbose)
        # No deflation (No residual matrices generated).
        Y <- NULL
        for (b in 1:J) Y[[b]] <- result$Y[,b, drop = FALSE]
        #Average Variance Explained (AVE) per block
        for (j in 1:J) AVE_X[[j]] =  mean(cor(A[[j]], Y[[j]])^2)
        
        #AVE outer 
        AVE_outer <- sum(pjs * unlist(AVE_X))/sum(pjs)
        
        AVE <- list(AVE_X = AVE_X, 
                    AVE_outer = AVE_outer,
                    AVE_inner = result$AVE_inner)
        
        a <- lapply(result$a, cbind)
        
        for (b in 1:J) {
          rownames(a[[b]]) = colnames(A[[b]])  
          rownames(Y[[b]]) = rownames(A[[b]])
          colnames(Y[[b]]) = "comp1"
        }
        
        out <- list(Y=Y, a=a, astar=a, 
                    C=C, scheme=scheme, c1=c1, ncomp=ncomp, 
                    crit = result$crit[length(result$crit)],
                    AVE = AVE)
        class(out) <- "sgcca"
        return(out)
    }
    
    ##################
    # Initialization #
    ##################
    
    Y <- NULL
    R <- A
    P <- a <- astar <- NULL
    crit <- list()
    AVE_inner <- rep(NA,max(ncomp))
    
    for (b in 1:J) P[[b]] <- a[[b]] <- astar[[b]] <- matrix(NA,pjs[[b]],N+1)
    for (b in 1:J) Y[[b]] <- matrix(NA,nb_ind,N+1)
    
    ##############################################
    #               If any ncomp > 1             #
    #      Determination of SGCCA components     #
    ##############################################
    
    
    
    for (n in 1:N) {
      if (verbose) cat(paste0("Computation of the SGCCA block components #", n, " is under progress... \n"))
      if(is.vector(c1)){
        sgcca.result <- sgccak(R, C, c1 = c1 , scheme=scheme, init = init, bias = bias, tol = tol, verbose=verbose) 
      } else{
        sgcca.result <- sgccak(R, C, c1 = c1[n, ] , scheme=scheme, init = init, bias = bias, tol = tol, verbose=verbose)
      }
      AVE_inner[n] <- sgcca.result$AVE_inner
      crit[[n]] <- sgcca.result$crit
      
      
      for (b in 1:J) Y[[b]][,n] <- sgcca.result$Y[ ,b]
      for (q in which(n <ndefl)) if(sum(sgcca.result$a[[q]]!=0) <= 1) warning(sprintf("Deflation failed because only one variable was selected for block #",q,"! \n"))
	    defla.result <- defl.select(sgcca.result$Y, R, ndefl, n, nbloc = J)
      R <- defla.result$resdefl
      for (b in 1:J) {
        P[[b]][,n] <- defla.result$pdefl[[b]]
        a[[b]][,n] <- sgcca.result$a[[b]]
      }
      
      if (n==1) {
        for (b in 1:J) astar[[b]][,n] <- sgcca.result$a[[b]]
      } 
      else {
        for (b in 1:J) astar[[b]][,n] <- sgcca.result$a[[b]] -  astar[[b]][,(1:n-1),drop=F] %*% drop( t(a[[b]][,n]) %*% P[[b]][,1:(n-1),drop=F] ) 
      }     
    }
    if (verbose) cat(paste0("Computation of the SGCCA block components #", N+1, " is under progress...\n"))
    if(is.vector(c1)) {
      sgcca.result <- sgccak(R, C, c1 = c1, scheme=scheme, init = init, bias = bias, tol = tol, verbose=verbose) 
    } else{
      sgcca.result <- sgccak(R, C, c1 = c1[N+1, ], scheme=scheme, init = init, bias = bias, tol = tol, verbose=verbose)
    }
    AVE_inner[max(ncomp)] <- sgcca.result$AVE_inner
    
    crit[[N+1]] <- sgcca.result$crit
    for (b in 1:J) {
      Y[[b]][,N+1]     <- sgcca.result$Y[, b]
      a[[b]][,N+1]     <- sgcca.result$a[[b]]
      astar[[b]][,N+1] <- sgcca.result$a[[b]] -  astar[[b]][,(1:N),drop=F] %*% drop( t(a[[b]][,(N+1)]) %*% P[[b]][,1:(N),drop=F] ) 
      rownames(a[[b]]) = rownames(astar[[b]]) = colnames(A[[b]])  
      rownames(Y[[b]]) = rownames(A[[b]])
      colnames(Y[[b]]) = paste0("comp", 1:max(ncomp))
    }
    
    shave.matlist <- function(mat_list, nb_cols) mapply(function(m, nbcomp) m[, 1:nbcomp, drop = FALSE], mat_list, nb_cols, SIMPLIFY=FALSE)
    shave.veclist <- function(vec_list, nb_elts) mapply(function(m, nbcomp) m[1:nbcomp], vec_list, nb_elts, SIMPLIFY=FALSE)    
    
    #Average Variance Explained (AVE) per block
    for (j in 1:J) AVE_X[[j]] =  apply(cor(A[[j]], Y[[j]])^2, 2, mean)
    
    #AVE outer 
    outer = matrix(unlist(AVE_X), nrow = max(ncomp))
    for (j in 1:max(ncomp)) AVE_outer[j] <- sum(pjs * outer[j, ])/sum(pjs)
    
    Y = shave.matlist(Y, ncomp)
    AVE_X = shave.veclist(AVE_X, ncomp)
    
    AVE <- list(AVE_X = AVE_X, 
                AVE_outer = AVE_outer,
                AVE_inner = AVE_inner)
    
    out <- list(Y = shave.matlist(Y, ncomp), 
                a = shave.matlist(a, ncomp), 
                astar = shave.matlist(astar, ncomp),
                C = C, c1 = c1, scheme = scheme,
                ncomp = ncomp, crit = crit,
                AVE = AVE
                )
    class(out) <- "sgcca"
    return(out)
    
}
