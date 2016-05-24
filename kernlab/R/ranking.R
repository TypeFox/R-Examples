## manifold ranking
## author: alexandros

setGeneric("ranking",function(x, ...) standardGeneric("ranking"))
setMethod("ranking",signature(x="matrix"),
          function (x,
                    y,
                    kernel    = "rbfdot",
                    kpar      = list(sigma = 1),
                    scale     = FALSE,
                    alpha     = 0.99,
                   iterations = 600,
                    edgegraph = FALSE,
                    convergence = FALSE,
                    ...)
          {
            m <- dim(x)[1]
            d <- dim(x)[2]
            if(length(y) != m)
              {
                ym <- matrix(0,m,1)
                ym[y] <- 1
                y <- ym
              }
            if (is.null(y))
              y <- matrix(1, m, 1) 
            labelled <- y != 0
            if (!any(labelled)) stop("no labels sublied")
            
            if(is.character(kernel))
              kernel <- match.arg(kernel,c("rbfdot","polydot","tanhdot","vanilladot","besseldot","laplacedot"))
            
           
            if(!is(kernel,"kernel"))
              {
                if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
                kernel <- do.call(kernel, kpar)
              }
    
            if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")

            if(scale)
             x <-  scale(x)
            ## scaling from ksvm     
            ## normalize ?
            
            if (is(kernel)[1]=='rbfkernel' && edgegraph){
              sigma = kpar(kernel)$sigma
              n <- dim(x)[1]
              dota <- rowSums(x*x)/2
              sed <- crossprod(t(x))
              for (i in 1:n)
                sed[i,] <-  - 2*(sed[i,] - dota - rep(dota[i],n))
              diag(sed) <- 0
              K <- exp(- sigma * sed)
             
              mst <- minimum.spanning.tree(sed)
              algo.mst <- mst$E
              max.squared.edge.length <-  mst$max.sed.in.tree
              edgegraph <- (sed <= max.squared.edge.length)
              K[!edgegraph] <- 0
              ##algo.edge.graph <- sparse(algo.edge.graph)
              rm(sed)
              gc()
            }
            else
              {
                edgegraph <- matrix()
                K <- kernelMatrix(kernel,x)   
              }
                        
            if (edgegraph && is(kernel)[1]!="rbfkernel"){
              warning('edge graph is only implemented for use with the RBF kernel')
              edgegraph <- matrix()
            }
            
            diag(K) <- 0
            ##K <- sparse(K)
            cs <- colSums(K)
            ##cs[cs <= 10e-6] <- 1
            
            D <- 1/sqrt(cs)
            K <- D * K %*% diag(D)
            
            if(sum(labelled)==1)
              y <- K[, labelled,drop = FALSE]
            else
              y <- as.matrix(colSums(K[, labelled]))
            K <- alpha * K[, !labelled]
            ym <- matrix(0,m,iterations)
            ym[,1] <- y 
            for (iteration  in  2:iterations)
              ym[, iteration] <- ym[, 1] + K %*% ym[!labelled, iteration-1]

            ym[labelled,] <- NA
            r <- ym
            r[!labelled,] <- compute.ranks(-r[!labelled, ])
            if(convergence)
              convergence <- (r - rep(r[,dim(r)[2]],iterations))/(m-sum(labelled))
            else
              convergence <- matrix()
            res <- cbind(t(t(1:m)), ym[,iterations], r[,iterations])
            return(new("ranking", .Data=res, convergence = convergence, edgegraph = edgegraph))
          })


## kernelMatrix interface
setMethod("ranking",signature(x="kernelMatrix"),
          function (x,
                    y,
                    alpha     = 0.99,
                   iterations = 600,
                    convergence = FALSE,
                    ...)
          {
            m <- dim(x)[1]

            if(length(y) != m)
              {
                ym <- matrix(0,m,1)
                ym[y] <- 1
                y <- ym
              }
            if (is.null(y))
              y <- matrix(1, m, 1) 
            labelled <- y != 0
            if (!any(labelled)) stop("no labels sublied")
                        
            diag(x) <- 0
            ##K <- sparse(K)
            cs <- colSums(x)
            ##cs[cs <= 10e-6] <- 1
            
            D <- 1/sqrt(cs)
            x <- D * x %*% diag(D)
            
            if(sum(labelled)==1)
              y <- x[, labelled,drop = FALSE]
            else
              y <- as.matrix(colSums(x[, labelled]))
            x <- alpha * x[, !labelled]
            ym <- matrix(0,m,iterations)
            ym[,1] <- y 
            for (iteration  in  2:iterations)
              ym[, iteration] <- ym[, 1] + x %*% ym[!labelled, iteration-1]

            ym[labelled,] <- NA
            r <- ym
            r[!labelled,] <- compute.ranks(-r[!labelled, ])
            if(convergence)
              convergence <- (r - rep(r[,dim(r)[2]],iterations))/(m-sum(labelled))
            else
              convergence <- matrix()
            res <- cbind(t(t(1:m)), ym[,iterations], r[,iterations])
            return(new("ranking", .Data=res, convergence = convergence))
          })


## list interface
setMethod("ranking",signature(x="list"),
          function (x,
                    y,
                    kernel    = "stringdot",
                    kpar      = list(length = 4, lambda = 0.5),
                    alpha     = 0.99,
                    iterations = 600, convergence = FALSE, ...)
          {
            m <- length(x)

            if(length(y) != m)
              {
                ym <- matrix(0,m,1)
                ym[y] <- 1
                y <- ym
              }

            if (is.null(y))
              y <- matrix(1, m, 1) 
            labelled <- y != 0
            if (!any(labelled)) stop("no labels sublied")
                        
            if(is.character(kernel))
              kernel <- match.arg(kernel,c("rbfdot","polydot","tanhdot","vanilladot","besseldot","laplacedot"))
            
            if(!is(kernel,"kernel"))
              {
                if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
                kernel <- do.call(kernel, kpar)
              }
    
            if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")
                        
            edgegraph <- matrix()
            K <- kernelMatrix(kernel,x)   

            diag(K) <- 0
            ##K <- sparse(K)
            cs <- colSums(K)
            ##cs[cs <= 10e-6] <- 1
            
            D <- 1/sqrt(cs)
            K <- D * K %*% diag(D)
            
            if(sum(labelled)==1)
              y <- K[, labelled,drop = FALSE]
            else
              y <- as.matrix(colSums(K[, labelled]))
            K <- alpha * K[, !labelled]
            ym <- matrix(0,m,iterations)
            ym[,1] <- y 
            for (iteration  in  2:iterations)
              ym[, iteration] <- ym[, 1] + K %*% ym[!labelled, iteration-1]

            ym[labelled,] <- NA
            r <- ym
            r[!labelled,] <- compute.ranks(-r[!labelled, ])
            if(convergence)
              convergence <- (r - rep(r[,dim(r)[2]],iterations))/(m-sum(labelled))
            else
              convergence <- matrix()
            res <- cbind(t(t(1:m)), ym[,iterations], r[,iterations])
            return(new("ranking", .Data=res, convergence = convergence, edgegraph = NULL))
          })

minimum.spanning.tree <- function(sed)
  { 
    max.sed.in.tree <- 0
    E <- matrix(0,dim(sed)[1],dim(sed)[2])
    n <- dim(E)[1]
    C <- logical(n)
    cmp <- sed
    diag(cmp) <- NA
    ans <- min(cmp, na.rm = TRUE) 
    i <- which.min(cmp)
    j <- i%/%n + 1
    i <- i%%n +1 
    
    for (nC  in  1:n) {
      cmp <- sed
      cmp[C,] <- NA
      cmp[,!C] <- NA
      if(nC == 1)
       {
         ans <- 1
         i <- 1
       }
      else{
        ans <- min(cmp, na.rm=TRUE) 
        i <- which.min(cmp)}
      j <- i%/%n + 1
      i <- i%%n + 1
      E[i, j] <- nC
      E[j, i] <- nC
      C[i] <- TRUE
      max.sed.in.tree <- max(max.sed.in.tree, sed[i, j])
    } 
    ## E <- sparse(E)
    res <- list(E=E, max.sed.in.tree=max.sed.in.tree)
  }

compute.ranks <- function(am) {
 
  rm <- matrix(0,dim(am)[1],dim(am)[2])
  for (j in 1:dim(am)[2])
    {
      a <- am[, j]
      sort <- sort(a, index.return = TRUE)
      sorted <- sort$x
      r <- sort$ix
      r[r] <- 1:length(r)
      
      while(1)
        {
          if(sum(na.omit(diff(sorted) == 0)) == 0)
            break
          tied <- sorted[min(which(diff(sorted) == 0))]
          sorted[sorted==tied] <- NA
          r[a==tied] <- mean(r[a==tied])
        } 
      rm[, j] <- r
    } 
  return(rm)
}

setMethod("show","ranking",
          function(object)
          {  cat("Ranking object of class \"ranking\"","\n")
             cat("\n")
             show(object@.Data)
             cat("\n")
             if(!any(is.na(convergence(object))))
               cat("convergence matrix included.","\n")
             if(!any(is.na(edgegraph(object))))
               cat("edgegraph matrix included.","\n")
           })
