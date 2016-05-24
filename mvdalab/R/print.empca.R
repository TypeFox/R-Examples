print.empca <- function(x, ncomp = x$impute.ncomps, ...) {
  if(is.null(x$CV.Results)) {
    print(x$Imputed.DataFrames[[ncomp]])
  } else {
    CVs <- x$CV.Results
    names(CVs) <- paste("Ncomp Imputation = ", 1:ncol(CVs))
    CVs$pca.ncomp <- 1:nrow(CVs)
    df <- melt(CVs, id = "pca.ncomp")
    graph.cvs <- with(df, ggplot(df, aes_string(as.factor(pca.ncomp), value)) + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      geom_point(aes_string(col = variable), size = 5) + 
      geom_line(aes_string(col = variable, group = variable), lwd = 1.25) + 
      ylab("GVC MSEP Approximation") + 
      xlab("PCA on Imputed Dataframe\nwith this No. of PCs") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(legend.title = element_blank()) + 
      ggtitle("Assessment of EM PCA") + 
      theme(strip.text.x = element_text(size = 15, face = "bold")) + 
      theme(plot.title = element_text(size = 20)) + 
      theme(axis.title.x = element_text(size = 20)) +
      theme(axis.title.y = element_text(size = 20, angle = 90)) + 
      theme(axis.text.x = element_text(size = 15, angle = 0, vjust = 0.5, face = "bold")) + 
      theme(axis.text.y = element_text(size = 15, angle = 0, face = "bold")))
    print(graph.cvs)
    print(x$Imputed.DataFrames[[ncomp]])
  }
}