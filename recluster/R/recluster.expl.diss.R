recluster.expl.diss<-function (tree, dist, maxcl = NULL, mincl=NULL, expld=TRUE) 
{ 
    dist <- as.matrix(dist)
    nclust <- NULL
    result <- NULL
    res <- NULL
    tree <- reorder(tree)
    mat <- nodeHeights(tree)
    if(is.null(mincl)){mincl<-1}
    if(is.null(maxcl)){maxcl<-nrow(mat)-1}
    mat2 <- rbind(c(0,1),mat[order(mat[, 1], mat[, 2]), ])
    mat2 <- mat2[mat2[, 1] + mat2[, 2] != 0, ]
    mat2 <- mat2[!duplicated(round(mat2[, 1],5)), ]
    if (maxcl>nrow(mat2)){maxcl<-nrow(mat2)-1}
    mat2 <- mat2[1:(maxcl+1), ]
    mat2 <- rbind(c(0,1),mat2[(mincl+1):nrow(mat2), ])
    matrix <- matrix(data = NA, ncol = nrow(mat2), nrow = length(tree$tip.label))
    comp <- rownames(dist)
    for (cl in 2:nrow(mat2)) {
        res <- treeSlice(tree, mat2[cl, 1] - 1e-06, trivial = TRUE)
        sub <- length(res)
        nclust[cl-1] <- sub
        for (subtrees in 1:sub) {
            tip <- length(res[[subtrees]]$tip.label)
            for (tp in 1:tip) {
                pos <- match(res[[subtrees]]$tip.label[tp], comp)
                matrix[pos, cl-1] <- subtrees
            }
        }
    }
    cluster <- NULL
    if(expld){
        beta <- sum(dist)    
        for (loops in 1:(nrow(mat2) - 1)) {
            cluster[loops] <-recluster.expl(dist, as.numeric(matrix[, loops]))
            }
      }   
   rownames(matrix) <- comp
    result$matrix <- matrix[,-ncol(matrix)]
    result$expl.div <- cluster
    result$nclust <- nclust
    return(result)
}

