XresidualContrib <- function(object, ncomp = object$ncomp, obs1 = 1) 
{
  X.residual.matrix <- as.matrix(object$Xdata) - 
    (object$scores[ , 1:ncomp] %*% t(object$loadings[ , 1:ncomp]))
  SPEc <- X.residual.matrix^2 * sign(X.residual.matrix)
  df <- data.frame(SPEc[obs1, ], Variable = names(SPEc[obs1, ]))
  names(df)[1] <- "Contributions"; row.names(df) <- NULL
  df$ncomp <- ncomp
  print(with(df, ggplot(df, aes_string(reorder(Variable, -Contributions, mean), Contributions, ymin = 0, ymax = Contributions)) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_linerange(aes_string(col = Variable), lwd = 2) + 
    geom_hline(yintercept = 0) + 
    theme(legend.position = "none") + 
    xlab("Variable") + 
    ylab("Contribution") + 
    facet_wrap(~ncomp) + 
    ggtitle("Contribution Plot for Xresiduals") + 
    theme(plot.title = element_text(size = 20)) + 
    theme(axis.title.x = element_text(size = 20)) + 
    theme(strip.text.x = element_text(size = 20, colour = "black", face="bold")) + 
    theme(strip.text.y = element_text(size = 20, colour = "black", face="bold")) + 
    theme(axis.title.y = element_text(size = 20, angle = 90)) + 
    theme(axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, face = "bold")) + 
    theme(axis.text.y = element_text(size = 15, angle = 0, face = "bold"))))
  df
}