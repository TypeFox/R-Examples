acfplot <- function(object, parm = NULL) {
  if(is.null(parm)) {  
    parm <- sample(names(object$Xdata), 1)
  } else parm <- parm
  Acf <- acf(smc(object)$error[, parm], plot = F)$acf[-1, 1, 1]
  clim0 <- qnorm((1 + .95)/2)/sqrt(length(Acf))
  clim <- clim0 * sqrt(cumsum(c(1, 2 * Acf^2)))
  ACF.LB <- max(-clim)
  ACF.UB <- min(clim)
  df <- data.frame(Seq = 1:length(Acf), acf = Acf)
  print(with(df, ggplot(df, aes_string(x = Seq, y = acf, ymin = 0, ymax = acf)) + 
          theme_bw() + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
          geom_linerange(lwd = 2) + 
          geom_hline(yintercept = 0) + 
          geom_hline(yintercept = ACF.UB, col = "red") + 
          geom_hline(yintercept = ACF.LB, col = "red") + 
          ggtitle(paste("Auto-correlation function on residuals for", parm, "in sMC Error Matrix")) +           
          theme(legend.position = "none") + 
          theme(plot.title = element_text(size = 20)) + 
          theme(axis.title.x = element_text(size = 20)) +
          theme(axis.title.y = element_text(size = 20, angle = 90)) + 
          theme(axis.text.x = element_text(size = 15, angle = 0, vjust = 0.5, face = "bold")) + 
          theme(axis.text.y = element_text(size = 15, angle = 0, face = "bold"))))
}
