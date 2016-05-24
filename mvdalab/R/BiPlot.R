BiPlot <- function(object, diag.adj = c(0, 0), axis.scaling = 2, cov.scale = FALSE) {
  comps <-  c(1, 2)
  if(object$method == "PCA") {
    D.. <- sqrt(diag(object$D) * (nrow(object$Xdata) - 1))[c(comps[1], comps[2])]
    U.. <- object$scores[, c(comps[1], comps[2])] %*% solve(diag(D..))
    V.. <- object$loadings[, c(comps[1], comps[2])]
  } else {
    D.. <- sqrt(diag(object$D2) * (nrow(object$Xdata) - 1))[c(comps[1], comps[2])]
    U.. <- object$scores[, c(comps[1], comps[2])] %*% solve(diag(D..))
    V.. <- object$loadings[, c(comps[1], comps[2])]
  }
  if(cov.scale == FALSE){
    G <- data.frame(U.. %*% diag(D..^diag.adj[1])); G$label<- 1:nrow(G)
    H <- data.frame(t(diag(D..^diag.adj[2]) %*% t(V..))); H$label<- row.names(H)
  } else {
    G <- data.frame(sweep(U.. %*% diag(D..^0), 2, sqrt(nrow(object$Xdata) - 1), "*")); G$label<- 1:nrow(G)
    H <- data.frame(sweep(t(diag(D..^0) %*% t(V..)), 2, sqrt(nrow(object$Xdata) - 1), "/")); H$label<- row.names(H)
  }
  names(G) <- names(H) <- c("A", "B", "label")
  par(mar=c(5.1, 7.1, 4.1, 7.1), xpd = F)
  plot(B ~ A, data = G, type = "n", 
       xlim = (c(min(G$B), max(G$B)) - mean(c(min(G$B), max(G$B)))) * axis.scaling, 
       ylim = (c(min(G$A), max(G$A)) - mean(c(min(G$A), max(G$A)))) * axis.scaling, col = "red", 
       xlab = "", ylab = "")
  title(main = "Biplot", line = 2.5, adj = 0)
  title(xlab = paste("PC ", comps[1], "Scores"), line = 2.5)
  title(ylab = paste("PC ", comps[2], "Scores"), line = 2.5)
  mtext(paste("PC ", comps[1], "Loadings"), side = 3, line = 2.5)
  mtext(paste("PC ", comps[2], "Loadings"), side = 4, line = 2.5)
  text(G$A, G$B, G$label, cex = 1, col = "red")
  par(new = TRUE)
  dev.hold()
  on.exit(dev.flush(), add = TRUE)
  c(min(H$B), max(H$B)) - mean(c(min(H$B), max(H$B)))
  plot(B ~ A, data = H, axes = FALSE, type = "n", 
       xlim = (c(min(H$B), max(H$B)) - mean(c(min(H$B), max(H$B)))) * axis.scaling, 
       ylim = (c(min(H$A), max(H$A)) - mean(c(min(H$A), max(H$A)))) * axis.scaling, col = "blue", 
       xlab = "", ylab = "")
  text(H$A, H$B, H$label, cex = 1, col = "blue")
  axis(3); axis(4)
  arrows(0, 0, H$A * 0.8, H$B * 0.8, col = "blue", lty = 2)
  par(mar = c(5, 4, 4, 2) + 0.1) 
}
    
    
    