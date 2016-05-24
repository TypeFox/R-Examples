PE <- function(object) {
  if(object$method == "bidiagpls" | object$method == "pls1gm") {
    XMSS <- diag(crossprod(object$scores %*% object$D2))
    XTSS <- sum(diag(crossprod(as.matrix(object$Xdata))))
    Ind <- data.frame(Percentages = (XMSS / XTSS), 
            Type = "Individual.Percent", Component = 1:object$ncomp)
    Cumulative <- data.frame(Percentages = cumsum(XMSS / XTSS), 
                  Type = "Cumulative.Percent", Component = 1:object$ncomp)
    df <- rbind(Ind, Cumulative)
    df$Type <- factor(df$Type, levels(df$Type)[c(2, 1)])
    print(with(df, ggplot(df, aes_string(x = Component, y = Percentages)) + 
            theme_bw() + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
            geom_line(aes_string(col = Type), lwd = 1.25) + 
            facet_wrap(~Type, scales = "free") + 
            ggtitle("Percent Explained of X") + 
            ylab("% Variance Explaned") + 
            xlab("Latent Variable Number") + 
            theme(legend.position = "none") + 
            theme(plot.title = element_text(size = 20)) + 
            theme(axis.title.x = element_text(size = 20)) +
            theme(strip.text.x = element_text(size = 20, colour = "black", face="bold")) + 
            theme(strip.text.y = element_text(size = 20, colour = "black", face="bold")) + 
            theme(axis.title.y = element_text(size = 20, angle = 90)) + 
            theme(axis.text.x = element_text(size = 15, angle = 0, vjust = 0.5, face = "bold")) + 
            theme(axis.text.y = element_text(size = 15, angle = 0, face = "bold")) + 
            scale_x_continuous(breaks = 1:object$ncomp)))
    return(df)
  } else {
df <- data.frame(stack(object$Percents.Explained[, 1:2]), 
                   object$Percents.Explained[, 3])
  names(df) <- c("Percentages", "Type", "Component")
  print(with(df, ggplot(df, aes_string(x = Component, y = Percentages)) + 
          theme_bw() + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
          geom_line(aes_string(col = Type), lwd = 1.25) + 
          facet_wrap(~Type, scales = "free") + 
          ggtitle("Percent Explained of X") + 
          ylab("% Variance Explaned") + 
          xlab("Latent Variable Number") + 
          theme(legend.position = "none") + 
          theme(plot.title = element_text(size = 20)) + 
          theme(axis.title.x = element_text(size = 20)) +
          theme(strip.text.x = element_text(size = 20, colour = "black", face="bold")) + 
          theme(strip.text.y = element_text(size = 20, colour = "black", face="bold")) + 
          theme(axis.title.y = element_text(size = 20, angle = 90)) + 
          theme(axis.text.x = element_text(size = 15, angle = 0, vjust = 0.5, face = "bold")) + 
          theme(axis.text.y = element_text(size = 15, angle = 0, face = "bold")) + 
          scale_x_continuous(breaks = object$Percents.Explained[, 3])))
  return(df)
  }
}
