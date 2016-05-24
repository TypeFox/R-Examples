plot.mvcomp <- function(x, Diff2Plot = c(3, 4), segments = 51, include.zero = FALSE, ...) {
  S2 <- x$Spooled[Diff2Plot, Diff2Plot]
  n1 <- x$n1
  n2 <- x$n2
  c2 <- x$c2
  center <- x$center
  Joint.CIs <- x$Joint.CIs[, -3]
  S2inv <- solve(S2)
  Radius <- sqrt(svd(S2)$d) * sqrt(((1 / n1) + (1 / n2 )) * c2)
  segments <- 51
  angles <- (0:segments) * 2 * pi/segments
  unit.circle <- cbind(cos(angles), sin(angles))
  Q <- sweep(svd(S2)$v, 2, Radius, "*")
  df <- data.frame(sweep(unit.circle %*% 
                               t(Q), 2, as.vector(center[Diff2Plot, ]), "+"))
  names(df) <- c("mu 1", "mu 2")
  output <- with(df, ggplot(df, aes_string(`mu 1`, `mu 2`)) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_path() + 
    theme(legend.position = "none") + 
    ggtitle("Confidence Ellipse for Mean Vector Differences") + 
    xlab(paste("mu", Diff2Plot[1], "difference")) + 
    ylab(paste("mu", Diff2Plot[2], "difference")) + 
    geom_vline(xintercept = Joint.CIs[Diff2Plot[1], 1], col = "lightgrey") + 
    geom_vline(xintercept = Joint.CIs[Diff2Plot[1], 2], col = "lightgrey") +
    geom_hline(yintercept = Joint.CIs[Diff2Plot[2], 1], col = "lightgrey") + 
    geom_hline(yintercept = Joint.CIs[Diff2Plot[2], 2], col = "lightgrey") +    
    theme(plot.title = element_text(size = 20)) + 
    theme(axis.title.x = element_text(size = 20)) + 
    theme(strip.text.x = element_text(size = 20, colour = "black", face="bold")) + 
    theme(strip.text.y = element_text(size = 20, colour = "black", face="bold")) + 
    theme(axis.title.y = element_text(size = 20, angle = 90)) + 
    theme(axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, face = "bold")) + 
    theme(axis.text.y = element_text(size = 15, angle = 0, face = "bold")))
  if(include.zero == TRUE) {
    output <- output + geom_vline(xintercept = 0, lty = 2) + geom_hline(yintercept = 0, lty = 2)
    print(output)
  } 
  print(output)
}
