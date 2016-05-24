coefsplot <- function(object, ncomp = object$ncomp, conf = .95) {
  if (object$val.method == "none" | object$val.method == "loo") {
    if(length(ncomp) == 1) {
      ncomp <- ncomp
      Stack <- data.frame(object$coefficients[, ncomp])
      names(Stack) <- ncomp
      df <- reshape(Stack, idvar = "variable", varying = list(1),
                      times = names(Stack[1]), 
                      v.names = "value", direction = "long", 
                      ids = row.names(Stack))
      row.names(df) <- NULL
      names(df)[1] <- "ncomp"
    } else {
      ncomp <- ncomp
      Stack <- data.frame(object$coefficients[, ncomp])
      names(Stack) <- ncomp
      Stack$variable <- row.names(Stack)
      df <- melt(Stack, measure = 1:length(ncomp))
      names(df)[2] <- "ncomp"
    }
    print(with(df, ggplot(df, aes(reorder(variable, -abs(value), mean), value)) + 
            theme_bw() + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
            geom_point(aes(col = variable), size = 5) + 
            geom_hline(yintercept = 0) + 
            ylab("coefficients") + 
            xlab("Variables") + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            ggtitle("Coefficients Plot") + 
            facet_wrap(~ncomp) + 
            theme(strip.text.x = element_text(size = 15, face = "bold")) + 
            xlab("Variables") + 
            theme(legend.position = "none") + 
            theme(plot.title = element_text(size = 20)) + 
            theme(axis.title.x = element_text(size = 20)) +
            theme(axis.title.y = element_text(size = 20, angle = 90)) + 
            theme(axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, face = "bold")) + 
            theme(axis.text.y = element_text(size = 15, angle = 0, face = "bold"))))
  } else {
    A <- do.call("rbind", coefficients.mvdareg(object, ncomp, conf = conf))
    df <- melt(A[, c(1, 6, 4:5, 2)], measure = 3:4, variable.name = "CI")
    df <- subset(df, ncomp = ncomp)
    print(with(df, ggplot(df, aes(reorder(variable, -abs(value), mean), value)) + 
            theme_bw() + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
            geom_line(aes(col = variable)) + 
            geom_hline(yintercept = 0) + 
            ylab("coefficients") + 
            xlab("Variables") + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            ggtitle("Coefficients Plot") + 
            facet_wrap(~ncomp) + 
            theme(strip.text.x = element_text(size = 15, face = "bold")) + 
            theme(legend.position = "none") + 
            theme(plot.title = element_text(size = 20)) + 
            theme(axis.title.x = element_text(size = 20)) +
            theme(axis.title.y = element_text(size = 20, angle = 90)) + 
            theme(axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, face = "bold")) + 
            theme(axis.text.y = element_text(size = 15, angle = 0, face = "bold"))))
    A
  }
}

