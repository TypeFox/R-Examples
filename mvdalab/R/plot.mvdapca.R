plot.mvdapca <- function(x, ...) {
  object <- x
  df <- data.frame(ncomp = 1:length(object$GVC), GVC = object$GVC)
  print(with(df, ggplot(df, aes_string(x = ncomp, y = GVC)) + 
               theme_bw() + 
               theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
               geom_line(lwd = 1.25) + 
               ggtitle("GVC Approximation of MSEP") + 
               xlab("Ncomps") + 
               ylab("MSEP") + 
               theme(legend.position = "none") + 
               theme(plot.title = element_text(size = 20)) + 
               theme(axis.title.x = element_text(size = 20)) +
               theme(strip.text.x = element_text(size = 20, colour = "black", face="bold")) + 
               theme(strip.text.y = element_text(size = 20, colour = "black", face="bold")) + 
               theme(axis.title.y = element_text(size = 20, angle = 90)) + 
               theme(axis.text.x = element_text(size = 15, angle = 0, vjust = 0.5, face = "bold")) + 
               theme(axis.text.y = element_text(size = 15, angle = 0, face = "bold")) + 
               scale_x_continuous(breaks = df$ncomp)))
}