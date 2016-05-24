plot.cp <- function(x, ncomp = "Overall", ...) {
  V <- x[[1]]
  df <- melt(cbind(subset(V, select = (colnames(V) %in% ncomp)), Variable = V$Variable), 
             id.vars = "Variable", value.name = "Contributions")
  print(with(df, ggplot(df, aes_string(reorder(Variable, -Contributions, mean), Contributions, ymin = 0, ymax = Contributions)) + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      geom_linerange(aes_string(col = Variable), lwd = 2) + 
      geom_hline(yintercept = 0) + 
      theme(legend.position = "none") + 
      xlab("Variable") + 
      ylab("Contribution") + 
      facet_wrap(~variable) + 
      ggtitle("Contribution Plot") + 
      theme(plot.title = element_text(size = 20)) + 
      theme(axis.title.x = element_text(size = 20)) + 
      theme(strip.text.x = element_text(size = 20, colour = "black", face="bold")) + 
      theme(strip.text.y = element_text(size = 20, colour = "black", face="bold")) + 
      theme(axis.title.y = element_text(size = 20, angle = 90)) + 
      theme(axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, face = "bold")) + 
      theme(axis.text.y = element_text(size = 15, angle = 0, face = "bold"))))
}