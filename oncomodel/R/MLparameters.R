`MLparameters` <-
function (x, tree, freq = NULL) 
{
    if (is.null(freq)) {
        naber <- ncol(x)
        dat <- cbind(x, x %*% 2^(0:(naber - 1)))
        x <- dat[which(!duplicated(dat[, naber + 1])), 1:naber]
        uniq <- unique(dat[, naber + 1])
        freq <- rep(NA, nrow(x))
        for (j in 1:length(freq)) freq[j] <- length(which(dat[, 
            naber + 1] == uniq[j]))
        
    }
    z <- 1
    bigtree <- tree
    currentsel <- 1:ncol(tree)
    parents <- match(tree[1, ], tree[2, ])
    siblings <- list()
    for (j in 1:ncol(tree)) siblings[[j]] <- which(tree[1, ] == 
        tree[1, j] & tree[2, ] != tree[2, j])
    repeat {
        z <- z + 1
        succ <- lapply(1:ncol(tree), function(x) which(tree[1, ] == tree[2, x]))
        indp <- indq <- matrix(0, nrow = nrow(x), ncol = ncol(tree))
        for (j in 1:nrow(x)) {
            activenodes <- which(x[j, ] == 1)
            while (any(activenodes %in% tree[2, ])) {
                edges <- match(activenodes, tree[2, ])
                edges <- edges[!is.na(edges)]
                indp[j, edges] <- 1
                activenodes <- unique(tree[1, edges])
            }
            for (k in 1:ncol(tree)) {
                if (indp[j, k] == 0 & (is.na(parents[k]) | indp[j, 
                  parents[k]] == 1)) 
                  indq[j, k] <- 1
            }
        }
        qvec <- freq %*% indq
        pvec <- freq %*% indp
        q2 <- rep(NA, ncol(tree))
        for (j in unique(tree[1, ])) {
            edges <- which(tree[1, ] == j)
            if (!any(tree[2, ] == j)) 
                q2[edges] <- qvec[edges]/(sum(freq))
            else {
                if (length(edges) == 2) {
                  q2[edges] <- qvec[edges]/pvec[edges[c(2, 1)]]
                }
                else {
                  loglik <- function(q) {
                    y <- pvec[edges] %*% log(1 - q) + qvec[edges] %*% 
                      log(q) - (pvec[edges[1]] + qvec[edges[1]]) * 
                      log(1 - prod(q))
                    return(-y)
                  }
                  result <- optim(par = rep(0.5, length(edges)), 
                    loglik, method = "L-BFGS-B", lower = rep(0.001, 
                      length(edges)), upper = rep(0.999, length(edges)))
                  q2[edges] <- result$par
                }
            }
        }
        p2 <- 1 - q2
        for (j in 1:ncol(tree)) {
            if (length(succ[[j]]) > 0) 
                p2[j] <- p2[j]/(1 - prod(q2[succ[[j]]]))
        }
        del <- which(p2 > 1.001)
        if (length(del) == 0) {
            pos <- which(q2 > 0)
            totloglik <- pvec[pos] %*% log(p2[pos]) + qvec[pos] %*% 
                log(q2[pos])
            break
        }
        else {
            currentsel <- currentsel[-del]
            for (k in seq(length(del), 1)) {
                if (!is.na(parents[del[k]])) 
                  tree[2, parents[del[k]]] <- tree[2, del[k]]
                if (length(siblings[[del[k]]]) > 0) 
                  tree[1, siblings[[del[k]]]] <- tree[2, del[k]]
                tree[which(tree >= tree[1, del[k]])] <- tree[which(tree >= 
                  tree[1, del[k]])] - 1
                tree <- tree[, -del[k]]
                parents <- match(tree[1, ], tree[2, ])
                siblings <- list()
                for (j in 1:ncol(tree)) siblings[[j]] <- which(tree[1, 
                  ] == tree[1, j] & tree[2, ] != tree[2, j])
            }
        }
    }
    smalltree <- tree
    tree <- bigtree
    p <- rep(1, ncol(tree))
    p[currentsel] <- p2
    output <- list(p, totloglik)
    names(output) <- c("p", "totloglik")
    return(output)
}

