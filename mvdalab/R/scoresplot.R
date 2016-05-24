scoresplot <- function(object, comps = c(1, 2), alphas = c(.95, .99), segments = 51) {
  if(ncol(object$scores) == 1) {
    df <- data.frame(A = object$scores[, comps[1]], label = 1:length(object$scores[, comps[1]]))
    print(with(df, ggplot(df, aes(x = label, y = A)) + 
          theme_bw() + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
          geom_point(aes(label = label, col = "black")) + 
          geom_line() + 
          ggtitle("Score Plot") + 
          theme(legend.position = "none") + 
          xlab(paste("Prin", comps[1])) + 
          ylim(-1, 1) + 
          geom_hline(yintercept = sd(df$A) + 2 * sd(df$A), col = "blue") + 
          geom_hline(yintercept = sd(df$A) - 2 * sd(df$A), col = "blue") + 
          geom_hline(yintercept = sd(df$A) + 3 * sd(df$A), col = "red") + 
          geom_hline(yintercept = sd(df$A) - 3 * sd(df$A), col = "red") + 
          ylab(paste("PC", comps[1])) + 
          theme(plot.title = element_text(size = 20)) + 
          theme(axis.title.x = element_text(size = 20)) +
          theme(axis.title.y = element_text(size = 20, angle = 90)) + 
          theme(axis.text.x = element_text(size = 15, angle = 0, vjust = 0.5, face = "bold")) + 
          theme(axis.text.y = element_text(size = 15, angle = 0, face = "bold"))))
    
  } else {
    df <- data.frame(A = object$scores[, comps[1]], B = object$scores[, comps[2]])
    df$label <- as.factor(1:nrow(object$Xdata))
    segments <- segments
    level <- alphas
    shape <- cor(df[, 1:2])
    center <- c(mean(df[, 1]),mean(df[, 2]))
    radius <- 1
    t. <- sapply(level, function(x) sqrt(qchisq(x, 2)))
    angles <- (0:segments) * 2 * pi/segments
    unit.circle <- cbind(cos(angles), sin(angles))
    Q <- diag(c(sd(df[, 1]),sd(df[, 2])))
    Ellipse <- llply(1:2, function(x) {
      Out <- data.frame(t(center + radius * t(t.[x] * unit.circle %*% Q)))
      names(Out) <- c("x", "y")
      Out
    })
    output.graph <- with(df, ggplot(df, aes(x = A, y = B)) + 
           theme_bw() + 
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
           geom_point() + 
           geom_text(aes(label = label, col = label)) + 
           ggtitle("Score Plot") + 
           theme(legend.position = "none") + 
           ylab(paste("PC", comps[2])) + 
           xlab(paste("PC", comps[1])) + 
           geom_hline(yintercept = 0) + 
           geom_vline(xintercept = 0) + 
           geom_path(data = Ellipse[[1]], aes(x = x, y = y), size = 1, linetype = 2, col = "blue") +
           geom_path(data = Ellipse[[2]], aes(x = x, y = y), size = 1, linetype = 2, col = "red") +       
           theme(plot.title = element_text(size = 20)) + 
           theme(axis.title.x = element_text(size = 20)) +
           theme(axis.title.y = element_text(size = 20, angle = 90)) + 
           theme(axis.text.x = element_text(size = 15, angle = 0, vjust = 0.5, face = "bold")) + 
           theme(axis.text.y = element_text(size = 15, angle = 0, face = "bold")))
 print(output.graph)
    }
  return(df)
}


