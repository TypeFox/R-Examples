weightsplot2D <- function(object, comps = c(1, 2)) {
  if(length(comps) == 1) {
    df <- data.frame(A = object$weights[, comps[1]], label = 1:length(object$weights[, comps[1]]))
    df$label <- row.names(df)
    row.names(df) <- NULL
    print(with(df, ggplot(df, aes(x = label, y = A, group = 1)) + 
            theme_bw() + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
            geom_point(aes(label = label, col = label, size = 5)) + 
            geom_line() + 
            ggtitle("Weights Plot") + 
            theme(legend.position = "none") + 
            xlab(paste("Prin", comps[1])) + 
            ylim(-1, 1) + 
            ylab("Weights") + 
            theme(plot.title = element_text(size = 20)) + 
            theme(axis.title.x = element_text(size = 20)) +
            theme(axis.title.y = element_text(size = 20, angle = 90)) + 
            theme(axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, face = "bold")) + 
            theme(axis.text.y = element_text(size = 15, angle = 0, face = "bold"))))
  } else {
    df <- data.frame(A = object$weights[, comps[1]], B = object$weights[, comps[2]])
    df$label <- row.names(df)
    row.names(df) <- NULL
    print(with(df, ggplot(df, aes(x = A, y = B, label = label, col = label)) + 
            theme_bw() + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
            geom_text(aes(label = label, col = label)) + 
            ggtitle("Weights Plot") + 
            theme(legend.position = "none") + 
            ylab(paste("Prin", comps[2])) + 
            xlab(paste("Prin", comps[1])) + 
            geom_segment(aes(x = 0, y = 0, xend = A, yend = B), 
                         arrow = arrow(), lwd = 1, lty = 2) + 
            geom_hline(yintercept = 0) + 
            geom_vline(xintercept = 0) + 
            xlab(paste("PC", comps[1])) + 
            ylab(paste("PC", comps[2])) + 
            theme(plot.title = element_text(size = 20)) + 
            theme(axis.title.x = element_text(size = 20)) +
            theme(axis.title.y = element_text(size = 20, angle = 90)) + 
            theme(axis.text.x = element_text(size = 15, angle = 0, vjust = 0.5, face = "bold")) + 
            theme(axis.text.y = element_text(size = 15, angle = 0, face = "bold"))))
  }
  return(df)
}
