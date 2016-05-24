bootstrap.size <-
function(bootstrap.grow.prune.results, plot.Ga=TRUE, filename=NULL, horizontal=TRUE){
  # EXTRACT NEEDED COMPONENTS
  boot.prune <- bootstrap.grow.prune.results$boot.prune
  tree0 <- bootstrap.grow.prune.results$initial.tree
  n <- tree0$size[1]
  
  OUT <- as.list(NULL)
  #  COMPUTE THE ALPHA PRIME'S
  prune0 <- boot.prune[[1]] 
  n.subtree <- nrow(prune0)
  alpha <- as.numeric(as.vector(prune0$alpha))
  # temp <- c(alpha[1], alpha[-length(alpha)])
  temp <- c(0, alpha[-length(alpha)])     #Set first value of alpha to 0
  alpha.prime <- sqrt(alpha*temp)
  # cbind(alpha,  alpha.prime=prune0$alpha.prime)
  b <- length(boot.prune)
  G <- as.numeric(as.vector(prune0$G))
  size.tmnl <- as.numeric(as.vector(prune0$size.tmnl))
  subtree <- as.numeric(as.vector(prune0$subtree))
  # tree.penalty <- log(nrow(teeth))
  G.a <- matrix(0, n.subtree, 5)
  penalty <- c(0, 2:4, log(n))
  for (i in 1:n.subtree) {
    a.i <- alpha.prime[i]
    bias <- 0
    for (j in 2:b){
      prune.bs <- boot.prune[[j]]
      alpha.bs <- as.numeric(as.vector(prune.bs$alpha))
      g <- as.numeric(as.vector(prune.bs$G))
      g.test <- as.numeric(as.vector(prune.bs$G.test))
      indx <- 1
      if (sum(alpha.bs <= a.i)>0) {
        temp1 <- which.max(which(alpha.bs<=a.i))
        indx <- ifelse(is.null(temp1), 1, temp1)
      }
      temp2 <- (g-g.test)[indx]
      bias <- bias + temp2
      # print(cbind(i, a.i, j, bias, indx, temp2))
    }
    G.honest <- G[i] - bias/(b-1) 
    G.a[i,] <- G.honest - penalty*(size.tmnl[i]-1)
  }
  out <- data.frame(cbind(size.tmnl, G.a))
  colnames(out) <- c("tmnl.size", "Ga", "Ga.2", "Ga.3", "Ga.4", "Ga.log_n")
  G.a <- out
  
  n.subtrees <- nrow(G.a)
  subtree.size <- G.a[,1]   
  # Plot Goodness of fit vs. a
  if (plot.Ga) {
    if (!is.null(filename)) postscript(file=filename, horizontal=horizontal)
    par(mfrow=c(1, 1), mar=rep(4, 4))
    
    x.min <- min(subtree.size); x.max <- max(subtree.size)
    y.min <- min(G.a$Ga.log_n); y.max <- max(G.a$Ga.2)
    plot(x=c(x.min, x.max), y=c(y.min, y.max), type="n", xlab="# of Terminal Nodes", ylab="G.a Values")
    for (j in 3:6) lines(subtree.size, G.a[,j], lty=1, col=j+4, lwd=2)
    legend(x.min+4, (y.min+y.max)/2, col=7:10, lty=1, lwd=2, legend=c("Ga.2", "Ga.3", "Ga.4", "Ga.ln(n)"))
    if (!is.null(filename)) dev.off()
  }
  #Obtain the best size tree
  bsize <- btree <- as.list(NULL)
  Ga.cols <- c("Ga.2", "Ga.3", "Ga.4", "Ga.log_n")
  for (j in Ga.cols) {
    best.size <- subtree.size[which.max(G.a[,j])]
    bsize[[j]] <- best.size
    btree[[j]] <- obtain.btree(tree0, bsize=best.size)
  }
  OUT$initial.tree <- tree0; OUT$G.a <- G.a; OUT$bsize <- bsize; OUT$btree <- btree
  return(OUT)
}
