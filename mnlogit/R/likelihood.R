#############################################################################
#     Loglikelihood function, gradient & Hessian evaluator.
#
# Args:
#   reponse - vector giving the response (for training data)
#   X - Design Matrix for Individual Specific Variable
#   Y - Design for choice attributes with individual variation & ch sp coeff
#   Z - Design for choice attributes and generic coeff
#   size - list with keys: N, K, p, f, d, nparams
#   coeffVec - Vector of coefficients (length = nparams)
#   Dimension specification for matrices:
#       See docs of newton.R
#   ncores  - number of processors allowed to use
#   hess - Evaluate Hessian & gradient only if TRUE
#   weights - a vector of frequency weights
#   loglikObj - If not NULL, gradient & the probability matrix are grabbed
#               from this, else they are computed
#
# Output:
#   loglik value, has attributes:
#      "probMat"  - probability Matrix
#      "gradient" - gradient vector
#      "hessian"  - hessian matrix (only is hess==TRUE) 
#      "funcTime" - Time spent in loglikelihood evaluation
#      "gradTime" - Time spent in gradient evaluation
#      "hessTime" - Time spent in Hessian evaluation
#     
# Notes:
#   Calls a C++ function to evaluate the Hessian.
#############################################################################
likelihood <- function(response, X, Y, Z, size, coeffVec, ncores, hess=TRUE, 
                       weights = NULL, loglikObj = NULL)
{
    t0 <- t1 <- t2 <- proc.time()[3]
    
    if (is.null(loglikObj)) {
        # Initialize utility matrix: dim(U) = N x K-1
        probMat <- matrix(rep(0, size$N * (size$K-1)), 
                          nrow = size$N, ncol = size$K-1)

        # First compute the utility matrix (stored in probMat)
        if (size$p) {
            probMat <- probMat + X %*% matrix(coeffVec[1:((size$K-1) *size$p)],
                                nrow = size$p, ncol = (size$K-1), byrow=FALSE)
        }
        if (size$f) {
            findYutil<- function(ch_k)
            {
                offset <- (size$K - 1)*size$p
                init <- (ch_k - 1)*size$N + 1
                fin <- ch_k * size$N
                Y[init:fin, , drop=FALSE] %*%
                  coeffVec[((ch_k-1)*size$f + 1 + offset):(ch_k*size$f+offset)]
            }
            vec <- as.vector(sapply(c(1:size$K), findYutil))
            vec <- vec - vec[1:size$N]
         
            probMat <- probMat + matrix(vec[(size$N+1):(size$N*size$K)], 
                                 nrow = size$N, ncol = (size$K-1), byrow=FALSE)
        }
        if (size$d) {
            probMat <- probMat +
              matrix(Z %*% coeffVec[(size$nparams - size$d + 1):size$nparams],
                     nrow = size$N, ncol=(size$K-1), byrow=FALSE)
        }

        # Compute partial log-likelihood
        loglik <- if (is.null(weights))
                     drop(as.vector(probMat) %*% response)
                  else 
                     drop(as.vector(probMat) %*% (weights * response))
        # Convert utility to probabilities - use logit formula
        probMat <- exp(probMat)                           # exp(utility)
        baseProbVec <- 1/(1 + rowSums(probMat))           # P_i0
        probMat <- probMat * matrix(rep(baseProbVec, size$K-1),
                          nrow = size$N, ncol = size$K-1) # P_ik
        # Negative log-likelihood
        loglik <- if (is.null(weights))
                      -1*(loglik + sum(log(baseProbVec)))
                  else
                      -1*(loglik + weights %*% log(baseProbVec))  

    } else {
        loglik <- loglikObj
        attributes(loglik) <- NULL
        probMat <- attr(loglikObj, "probMat")  
        baseProbVec <- 1 - rowSums(probMat)
    }
    t1 <- proc.time()[3]  # Time after function evaluation
    
    # Gradient calculation
    gradient <- if (hess) {
        responseMat <- matrix(response, nrow=size$N, ncol=(size$K-1))
        baseResp <- rep(1, size$N) - rowSums(responseMat)
   
        xgrad <- if (!is.null(X)) {
            if (is.null(weights))
                as.vector(crossprod(X, responseMat - probMat))
            else
                as.vector(crossprod(X, weights * (responseMat - probMat)))
        } else NULL
        ygrad <- if (!is.null(Y)) {
            yresp <- c((baseResp - baseProbVec), (response-as.vector(probMat)))
            findYgrad <- function(ch_k)
            {
              init <- (ch_k - 1)*size$N + 1
              fin <- ch_k * size$N
              if (is.null(weights))
                crossprod(Y[init:fin, , drop=FALSE], yresp[init:fin])
              else
                crossprod(Y[init:fin, , drop=FALSE], weights * yresp[init:fin])
            }
            as.vector(sapply(c(1:size$K), findYgrad))
        } else NULL
       
        if (!is.null(Z)) {
            zgrad <- if (is.null(weights))
                       as.vector(Z) * as.vector(responseMat - probMat)
                     else
                       as.vector(Z) * weights * as.vector(responseMat - probMat)
            zgrad <- colSums(matrix(zgrad, nrow=nrow(Z), ncol=size$d))
        } else zgrad <- NULL
        -1 * c(xgrad, ygrad, zgrad)
    }
    else  attr(loglikObj, "gradient")
    t2 <- proc.time()[3]  # Time after gradient evaluation

    # Hessian calculation
    ans <- if (hess) {
       hessMat <- rep(0, size$nparams * size$nparams)
       .Call("computeHessianDotCall" , as.integer(size$N), as.integer(size$K),
            as.integer(size$p), as.integer(size$f), as.integer(size$d), 
            if (is.null(X)) NULL else as.double(t(X)), 
            if (is.null(Y)) NULL else as.double(t(Y)), 
            if (is.null(Z)) NULL else as.double(Z),
            if (is.null(weights)) NULL else as.double(weights),
            probMat, baseProbVec, as.integer(ncores), hessMat)
       hessMat
    } else NULL
    t3 <- proc.time()[3]  # Time after Hessian evaluation
    
    # Prepare to return
    attr(loglik, "probMat")  <- probMat
    attr(loglik, "gradient") <- gradient
    attr(loglik, "hessian")  <- if (!hess) attr(loglikObj, "hessian") 
        else matrix(ans, nrow = size$nparams, ncol = size$nparams)
        #else matrix(ans$hessMat, nrow = size$nparams, ncol = size$nparams)
    attr(loglik, "funcTime") <- t1 - t0 
    attr(loglik, "gradTime") <- t2 - t1 
    attr(loglik, "hessTime") <- t3 - t2
    return(loglik) 
}
