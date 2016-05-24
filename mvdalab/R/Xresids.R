Xresids <- function(object, ncomp = object$ncomp, conf = c(.95, .99), normalized = TRUE) {
  X.residual.matrix <- as.matrix(object$Xdata) - 
  (object$scores[ , 1:ncomp] %*% t(object$loadings[ , 1:ncomp]))
  df.Score <- (apply(X.residual.matrix, 1, function(x) sum(x^2)))
  v <- var(df.Score)
  m <- mean(df.Score)
  Mult <- v / (2 * m)
  dfreed <- (2 * m^2) / v
  Upper.95 <- Mult * qchisq(conf[1], dfreed) 
  Upper.99 <- Mult * qchisq(conf[2], dfreed) 
  df <- data.frame(Seq = 1:length(df.Score), df.Score, Upper.95, Upper.99)
  df$Result.95 <- ifelse(df.Score > Upper.95 & df.Score < Upper.99, "Out", "in")
  df$Result.99 <- ifelse(df.Score > Upper.99, "Out", "in")
  df$Result.95 <- ifelse(df$Result.99 == "Out", "Out", df$Result.95)
  print(with(df, ggplot(df, aes_string(Seq, df.Score)) + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
        geom_text(aes_string(label = as.factor(Seq), col = as.factor(Seq))) +
        geom_line() + 
        geom_hline(yintercept = df$Upper.95, col = "blue") + 
        geom_hline(yintercept = df$Upper.99, col = "red") + 
        theme(legend.position = "none") + 
        xlab("Batch Sequence") + 
        ylab("Residuals") + 
        ggtitle(paste("X Residuals Range Plot for", ncomp, "components")) +           
        theme(plot.title = element_text(size = 20)) + 
        theme(axis.title.x = element_text(size = 20)) +
        theme(axis.title.y = element_text(size = 20, angle = 90)) + 
        theme(axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, face = "bold")) + 
        theme(axis.text.y = element_text(size = 15, angle = 0, face = "bold"))))
  return(df)
}