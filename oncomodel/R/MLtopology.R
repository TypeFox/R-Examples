`MLtopology` <-
function (x, verbose = FALSE) 
{
    if (!is.data.frame(x) & !is.matrix(x)) 
        stop("'x' should be of type 'matrix'.")
    x <- as.matrix(x)
    if (!(all(x %in% c(0, 1)))) 
        stop("'x' should be binary.")
    if (any(colSums(x) == 0)) 
        stop("Every event should occur in at least one case.")
    if (is.null(colnames(x)) | any(duplicated(colnames(x)))) 
        stop("'x' should have column names.")
    used <- sample(ncol(x), 2)
    oldtree <- cbind(c(3, 1), c(3, 2), c(4, 3))
    order <- sample((1:ncol(x))[-used])
    for (k in 1:length(order)) {
        oldtree[which(oldtree > length(used))] <- oldtree[which(oldtree > 
            length(used))] + 1
        bestvalue <- -Inf
        used <- c(used, order[k])
        dat <- x[, used]
        naber <- ncol(dat)
        dat <- cbind(dat, dat %*% 2^(0:(naber - 1)))
        newdata <- dat[which(!duplicated(dat[, naber + 1])), 
            1:naber]
        uniq <- unique(dat[, naber + 1])
        freq <- rep(NA, nrow(newdata))
        for (j in 1:length(freq)) freq[j] <- length(which(dat[, 
            naber + 1] == uniq[j]))
        tree <- cbind(c(2 * length(used), length(used)), c(2 * 
            length(used), oldtree[2, 1]), c(oldtree[1, 1], 2 * 
            length(used)), oldtree[, 2:ncol(oldtree)])
        res <- MLparameters(newdata, tree, freq)
        besttree <- tree
        if (!is.na(res$totloglik)) 
            bestvalue <- res$totloglik
        for (edge in 2:ncol(oldtree)) {
            tree <- cbind(oldtree[, which(1:ncol(oldtree) < edge)], 
                c(2 * length(used), length(used)), c(2 * length(used), 
                  oldtree[2, edge]), c(oldtree[1, edge], 2 * 
                  length(used)), oldtree[, which(1:ncol(oldtree) > 
                  edge)])
            if (!is.na(res$totloglik) & res$totloglik > bestvalue) {
                besttree <- tree
                bestvalue <- res$totloglik
            }
        }
        oldtree <- besttree
        if (verbose == TRUE) {
            cat("adding", colnames(x)[order[k]], "\t")
            cat(bestvalue, "\n")
        }
        oldbestvalue <- bestvalue
        curedge <- which(besttree[2, ] == length(used))
        repeat {
            flag <- 0
            for (rmedge in (1:ncol(oldtree))[-curedge]) {
                tree <- oldtree
                subtreeind <- c(subtree(tree[2, rmedge], tree), 
                  rmedge)
                if (ncol(tree) - length(subtreeind) > 2) {
                  subt <- as.matrix(tree[, subtreeind])
                  rest <- as.matrix(tree[, -subtreeind])
                  restparents <- match(rest[1, ], rest[2, ])
                  sib <- which(rest[1, ] == tree[1, rmedge])
                  if (length(sib) == 1) {
                    rest[1, sib] <- rest[1, restparents[sib]]
                    rest <- as.matrix(rest[, -restparents[sib]])
                  }
                  for (insedge in 1:ncol(rest)) {
                    insedgeintree <- which(oldtree[1, ] == rest[1, 
                      insedge] & oldtree[2, ] == rest[2, insedge])
                    subtreeind2 <- subtree(rest[2, insedge], 
                      rest)
                    subt2 <- as.matrix(rest[, subtreeind2])
                    rest2 <- rest[, -c(insedge, subtreeind2)]
                    newtree <- cbind(subt, subt2, c(oldtree[1, 
                      rmedge], rest[2, insedge]), c(rest[1, insedge], 
                      oldtree[1, rmedge]), rest2)
                    tree <- newtree
                    res <- MLparameters(newdata, tree, freq)
                    if (!is.na(res$totloglik) & res$totloglik - 
                      bestvalue > 1e-10) {
                      besttree <- tree
                      bestvalue <- res$totloglik
                      curedge <- ncol(subt)
                      flag <- 1
                    }
                  }
                }
            }
            if (flag == 1) {
                if (verbose == TRUE) 
                  cat("change due to rearrangement\t", bestvalue, 
                    "improvement by", bestvalue - oldbestvalue, 
                    "\n")
                oldtree <- besttree
                oldbestvalue <- bestvalue
            }
            else {
                break
            }
        }
    }
    totloglik <- bestvalue
    tree <- besttree
    res <- MLparameters(newdata, tree, freq)
    cat("resulting tree", res$totloglik, "\n")
    tr <- tree[2, ]
    tr[which(tree[2, ] <= naber)] <- used[tree[2, which(tree[2, 
        ] <= naber)]]
    tree[2, ] <- tr
    p <- res$p
    totloglik <- res$totloglik
    ##convert into Newick format:
    nodes <- sort(unique(as.vector(tree)))
    for (j in 1:length(nodes)) tree[which(tree == nodes[j])] <- j
    res <- as.list(c(colnames(x), rep("(", length(nodes) - ncol(x))))
    done <- rep(FALSE, length(nodes))
    done[1:ncol(x)] <- TRUE
    while (!all(done)) {
      waiting <- numeric()
      for (j in which(!done)) {
        if (all(done[tree[2, which(tree[1, ] == j)]])) 
          waiting <- c(waiting, j)
      }
      for (j in waiting) {
        edges <- which(tree[1, ] == j)
        for (k in 1:length(edges)) {
          tmp <- paste(res[[tree[2, edges[k]]]], ":", format(-log(p[edges[k]]), 
                                                             digits = 3), sep = "")
          res[[j]] <- paste(res[[j]], tmp, sep = "")
          if (k == length(edges)) 
            res[[j]] <- paste(res[[j]], ")", sep = "")
          else res[[j]] <- paste(res[[j]], ",", sep = "")
        }
        done[j] <- TRUE
      }
    }
    parents <- match(tree[1, ], tree[2, ])
    root <- tree[1, which(is.na(parents))]
    newick <- paste(res[[root]], ";", sep="")
    
    output <- list(tree, p, totloglik, colnames(x), newick)
    names(output) <- c("tree", "p", "totloglik", "var.names", "newick")
    return(output)
}

