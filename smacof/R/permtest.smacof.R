## MDS permutation test
permtest.smacof <- function(object, data,  method.dat = "pearson", nrep = 100, verbose = TRUE, ...) {

## val ... stress value  
## n... number of objects
## p... number of dimensions
## method ... "full" or "rows"  
    
    method.dat <- match.arg(method.dat, c("pearson", "spearman", "kendall", "euclidean", "maximum", "manhattan", "canberra", "binary"))
    #method.diss <- match.arg(method.diss, c("full", "rows"))
    method.diss <- "full"
    
    n <- object$nobj          ## number of objects
    if (!missing(data)) if(ncol(data) != n) stop("Number of columns need to match number of MDS objects!") 
  
    p <- object$ndim          ## number of dimensions
    val <- object$stress      ## stress-1 value
    dissvec <- as.vector(object$delta)  ## observed dissimilarities as vector
    smacall <- object$call
    
    m <- n*(n-1)/2            ## number of lower triangular elements
    irep <- 1
    str <- rep(0, nrep)      ## vector for stress values

    if (missing(data)) {                      ## permutations based on dissimilarit matrix
      repeat 
      {
        
        if (method.diss == "full") { 
          ## --- begin old permutation scheme
          delta <- matrix(0, n, n)
          delta[outer(1:n, 1:n, ">")] <- sample(dissvec, m)           ## sample dissimilarity matrix
          delta <- delta + t(delta)
          ## --- end old permutation scheme
        } else {
          ## --- alternative permutation scheme
          delta1 <- as.matrix(object$delta)
          permind <- sample(1:nrow(delta1))         ## index permutation row/columns
          delta <- delta1[permind, ]
          rownames(delta) <- colnames(delta) <- rownames(delta1)
        } 
        ## --- end alternative permutation scheme
        
        smacall$delta <- delta
        smRes <- eval(smacall)
        #smRes <- smacofSym(delta, ndim = p, ...)        
        str[irep] <- smRes$stress                               ## store stress of no-structure matrix 
        if (verbose) cat("Permutation: ", formatC (irep, digits=3, width=3), "Stress: ", formatC (str[irep], digits=10, width=15, format="f"), "\n")
        if (irep == nrep) break()
        irep <- irep + 1  
      }
    } else {                                ## permutations based on raw data
      
      N <- nrow(data)
      for (irep in 1:nrep) {
        proxmat <- diag(1, n, n)
         
         for (i in 2:n) {                                
           colperm <- data[sample(1:N),i:n]  
           if (any(method.dat == c("pearson", "spearman", "kendall"))) {
             proxmat[i-1, i:n] <- cor(data[,i-1], colperm, method = method.dat, use = "complete.obs")    ## compute proximities (correlation)
           } else {
             proxmat[i-1, i:n] <- as.matrix(dist(t(cbind(data[,i-1], colperm)), method = method.dat))[2:(n-i+2),1] ## compute dissimilarity
           }
         }

        proxmat[lower.tri(proxmat)] <- t(proxmat)[lower.tri(proxmat)] 
        
        if (any(method.dat == c("pearson", "spearman", "kendall"))) dissmat <- sim2diss(proxmat, ...) else dissmat <- proxmat
        smacall$delta <- dissmat
        smRes <- eval(smacall)  
        str[irep] <- smRes$stress                      ## MDS fit
        if (verbose) cat("Permutation: ", formatC (irep, digits=3, width=3), "Stress: ", formatC (str[irep], digits=10, width=15, format="f"), "\n")
        
      }
    }
    
    pval <- length(which(str < val))/nrep
      
    result <- list(stressvec = str, stress.obs = val, pval = pval, nobj = n, nrep = nrep, 
                   call = match.call())
    class(result) <- "smacofPerm"
    result
}