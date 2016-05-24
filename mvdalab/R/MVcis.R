MVcis <- function(data, segments = 51, level = .95, Vars2Plot = c(1, 2), 
                   include.zero = F) {
  X <- my.dummy.df(data)
  x <- as.matrix(colMeans(X), ncol = 1)
  n <- nrow(X)
  p <- ncol(X)
  S <- cov(X) / n
  Sinv <- solve(S)
  c2 <- sqrt((((n - 1) * p) / (n - p)) * qf(.95, p, n - p))
  Joint.CIs <- cbind(x - c2 * sqrt(apply(X, 2, function(x) var(x) / n)), 
                     x + c2 * sqrt(apply(X, 2, function(x) var(x)) / n))
  S2 <- S[Vars2Plot, Vars2Plot]
  S2inv <- solve(S2)
  Radius <- sqrt(svd(S2)$d) * c2
  segments = 51
  angles <- (0:segments) * 2 * pi/segments
  unit.circle <- cbind(cos(angles), sin(angles))
  Q <- sweep(svd(S2)$v, 2, Radius, "*")
  df <- data.frame(sweep(unit.circle %*% 
                                 t(Q), 2, as.vector(x[Vars2Plot, ]), "+"))
  names(df) <- c("mu 1", "mu 2")
  output <- with(df, ggplot(df, aes_string(`mu 1`, `mu 2`)) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_path() + 
    theme(legend.position = "none") + 
    ggtitle("Confidence Ellipse for Mean Vector Differences") + 
    xlab(rownames(x)[Vars2Plot[1]]) + 
    ylab(rownames(x)[Vars2Plot[2]]) + 
    geom_vline(xintercept = Joint.CIs[Vars2Plot[1], ], col = "lightgrey") + 
    geom_hline(yintercept = Joint.CIs[Vars2Plot[2], ], col = "lightgrey") +
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
  print(Joint.CIs)
}

