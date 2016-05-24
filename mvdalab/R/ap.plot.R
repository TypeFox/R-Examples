ap.plot <- function(object, ncomp = object$ncomp) {
  if(class(object) != "mvdareg") {
    stop("Coefficient only implemented for PLS")
  }
  ap.dat <- data.frame(Seq = 1:length(as.vector(object$Yactual)), 
                       'Actual' = as.vector(object$Yactual), 
                       Predicted = object$iPreds[, ncomp], 
                       'Residuals' = object$residuals[, ncomp])
  df <- melt(ap.dat, measure = c(2, 4), value.name = "values", variable.name = "ind")
  print(with(df, ggplot(df, aes_string(Predicted, values)) + 
          theme_bw() + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
          geom_text(aes_string(label = as.factor(Seq), col = as.factor(Seq))) + 
          geom_smooth(method = "lm", se = F) + 
          facet_wrap(~ind, scales = "free") + 
          xlab("Predicted") + 
          ylab("") + 
          theme(legend.position = "none") + 
          ggtitle(paste("Fit Diagnostic for", ncomp, "components")) +           
          theme(strip.text.x = element_text(size = 20, angle = 0)) + 
          theme(plot.title = element_text(size = 20)) + 
          theme(axis.title.x = element_text(size = 20)) +
          theme(axis.title.y = element_text(size = 20, angle = 90)) + 
          theme(axis.text.x = element_text(size = 15, angle = 0, vjust = 0.5, face = "bold")) + 
          theme(axis.text.y = element_text(size = 15, angle = 0, face = "bold"))))
  ap.dat
}



