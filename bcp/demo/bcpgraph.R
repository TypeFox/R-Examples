if(require("igraph")) {
  set.seed(1)
  
  # make the adjacency structure (for a short-boundary model)
  p <- 0.5
  adj <- matrix(0, 20, 20)
  adj[1:10,1:10] <- rbinom(100, 1, p)
  adj[11:20, 11:20] <- rbinom(100, 1, p)
  adj[5,15] <- adj[6,18] <- adj[2, 13] <- 1
  adj[lower.tri(adj)] <- t(adj)[lower.tri(adj)]
  diag(adj) <- 0
  g <- graph.adjacency(adj, mode='undirected')
  
  # generate the data
  z <- rep(c(0, 3), each=10)
  y <- z + rnorm(20)
  y.rounded <- round(y,1)
  
  # plot true means
  y.scale <- as.character(seq(-1,4, by=0.1))
  cols <- rev(heat.colors(length(y.scale)))
  lay <-layout.fruchterman.reingold(g, niter=100000)
  V(g)$name <- round(z,1)
  V(g)$color <- cols[match(as.character(round(z,1)), y.scale)]
  plot(g, edge.width=4, layout=lay, vertex.size=30, vertex.label.cex=1.5, main="True Means")
  lines(c(-1.3, 1.3), c(0.1,0.1), lwd=4, col="blue")
  
  # plot the data
  V(g)$name <- round(y,1)
  V(g)$color <- cols[match(as.character(round(z,1)), y.scale)]
  plot(g, edge.width=4, layout=lay, vertex.size=30, vertex.label.cex=1.5, main="Data")
  lines(c(-1.3, 1.3), c(0.1,0.1), lwd=4, col="blue")
  
  # format the adjacency information
  adj2 <- list()
  for (i in 1:20) {
    adj2[[i]] <- which(adj[i,] ==1) -1
  }
  
  # run bcp
  a <- bcp(y, p0=0.2, adj=adj2, burnin=1000, mcmc=500)
  
  ## plot of posterior means
  V(g)$name <- round(a$posterior.mean[,1],1)
  tmp.rounded <- round(a$posterior.mean[,1],1)
  V(g)$color <- cols[match(as.character(tmp.rounded), y.scale)]
  plot(g, edge.width=4, layout=lay, vertex.size=30, main="Posterior Means",
       vertex.label.cex=1.5)
  
  
  ## plot of posterior probability of lying at a block boundary
  V(g)$name <- round(a$posterior.prob,2)
  tmp.rounded <- round(a$posterior.prob,2)
  p.scale <- as.character(seq(0, 1, by=0.01))
  V(g)$color <- cols[match(as.character(tmp.rounded), p.scale)]
  plot(g, edge.width=4, layout=lay, vertex.size=30, main="Posterior Boundary Probabilities",
       vertex.label.cex=1.5)
}