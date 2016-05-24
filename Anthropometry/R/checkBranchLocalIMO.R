checkBranchLocalIMO <- function(tree,data,i,maxsplit,asw.tol,local.const,
                                orness,type,ah,verbose,...){
  if (is.vector(data)){
   proposal <- list(reject = TRUE, tree = -1)
  } else {
     if (ncol(data) <= 2){
      proposal <- list(reject = TRUE,tree = -1)
     } else {
        if(sum(tree$clustering == i) <= 2){ #First stopping criteria.
         proposal <- list(tree = tree,reject = TRUE)
        } else {
           which.x <- (tree$clustering == i)
           xi <- data[which.x,]

           if (nrow(xi) <= maxsplit){ 
            maxsplit2 <- max(nrow(xi) - 1, 2)
            DIST <- ext.dist(xi, maxsplit2, orness, ah, verbose)
            out <- INCAnumclu(DIST, K = maxsplit2, method = "pam", L = NULL, noise = NULL)
            maxsplit <- maxsplit2
           }else{
             DIST <- ext.dist(xi, maxsplit, orness, ah, verbose)
             out <- INCAnumclu(DIST, K = maxsplit, method = "pam", L = NULL, noise = NULL)
            }

           if(max(out$INCAindex[2:maxsplit]) <= 0.2){ #Second stopping criteria.
            proposal <- list(tree = tree,reject = TRUE) 
           }else{
             xi.ps <- getBestPamsamIMO(xi, maxsplit, orness = orness, type, 
                                       ah, verbose, ...)
                        
            if (is.null(local.const)){
             n.sub.clust <- xi.ps$num.of.clusters
             asw.vec <- rep(NA,n.sub.clust)
             
              for (j in 1:n.sub.clust){
               if (sum(xi.ps$clustering==j) <= 2){
                  asw.vec[j] <- 0
               } else {
                  xij <- xi[xi.ps$clustering==j,]
                  asw.vec[j] <- getBestPamsamIMO(xij, maxsplit, orness = orness, type, 
                                                 ah, verbose, ...)$asw
                 }
              }
                  if (xi.ps$asw > mean(asw.vec) - asw.tol){ #Third stopping criteria.
                   tree <- update.tree.local(object = tree, xi.ps, which.x, i) 
                   proposal <- list(tree = tree, reject = FALSE) 
                  } else {
                     proposal <- list(tree = tree, reject = TRUE)
                   }
             } else {
                if (xi.ps$asw > local.const){
                 tree <- update.tree.local(object = tree, xi.ps, which.x, i)
                 proposal <- list(tree = tree, reject = FALSE)
                } else {
                  proposal <- list(tree = tree, reject = TRUE)
                  }
               }
           }
         }
       }
     } 
    proposal
}








