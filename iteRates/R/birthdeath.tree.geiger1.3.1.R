birthdeath.tree.geiger1.3.1 <-
function (b, d, time.stop = 0, taxa.stop = 0, seed = 0, print.seed = FALSE, 
    return.all.extinct = TRUE) 
{
    if (seed == 0){
    	date<-date()
    	seed1 = as.numeric(strsplit(substring(date,12,19),":")[[1]])%*%c(1,100,10000)
    	seed <- runif(1, min=0, max=50) * seed1
    	set.seed(seed)
    	}
#        seed = set.seed.clock(print = print.seed)
    if (time.stop == 0 & taxa.stop == 0) 
        stop("Must have stopping criterion\n")
    while (1) {
        edge <- rbind(c(1, 2), c(1, 3))
        edge.length <- rep(NA, 2)
        stem.depth <- numeric(2)
        alive <- rep(TRUE, 2)
        t <- 0
        next.node <- 4
        repeat {
            if (taxa.stop) 
                if (sum(alive) >= taxa.stop) 
                  break
            if (sum(alive) == 0) 
                break
            dt <- rexp(1, sum(alive) * (b + d))
            t <- t + dt
            if (time.stop) 
                if (t >= time.stop) {
                  t <- time.stop
                  break
                }
            r <- runif(1)
            if (r <= b/(b + d)) {
                random_lineage <- round(runif(1, min = 1, max = sum(alive)))
                e <- matrix(edge[alive, ], ncol = 2)
                parent <- e[random_lineage, 2]
                alive[alive][random_lineage] <- FALSE
                edge <- rbind(edge, c(parent, next.node), c(parent, 
                  next.node + 1))
                next.node <- next.node + 2
                alive <- c(alive, TRUE, TRUE)
                stem.depth <- c(stem.depth, t, t)
                x <- which(edge[, 2] == parent)
                edge.length[x] <- t - stem.depth[x]
                edge.length <- c(edge.length, NA, NA)
            }
            else {
                random_lineage <- round(runif(1, min = 1, max = sum(alive)))
                edge.length[alive][random_lineage] <- t - stem.depth[alive][random_lineage]
                alive[alive][random_lineage] <- FALSE
            }
        }
        if (return.all.extinct == T | sum(alive) > 1) 
            break
    }
    edge.length[alive] <- t - stem.depth[alive]
    n <- -1
    for (i in 1:max(edge)) {
        if (any(edge[, 1] == i)) {
            edge[which(edge[, 1] == i), 1] <- n
            edge[which(edge[, 2] == i), 2] <- n
            n <- n - 1
        }
    }
    edge[edge > 0] <- 1:sum(edge > 0)
    tip.label <- 1:sum(edge > 0)
    mode(edge) <- "character"
    mode(tip.label) <- "character"
    obj <- list(edge = edge, edge.length = edge.length, tip.label = tip.label)
    class(obj) <- "phylo"
    obj <- old2new.phylo(obj)
    obj <- read.tree(text = write.tree(obj))
    obj
    #birthdeath.tree function from geiger v1.3-1
}

