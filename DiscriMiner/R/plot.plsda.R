#' @S3method plot plsda
plot.plsda <-
function(x, ...)
{
  ## Plot circle of correlations between variables
  # get correlations
  cor_tx = x$comp_vars
  cor_ty = x$comp_group
  # points for generating circle
  z = seq(0, 2*pi, l=100)
  # open plot
  plot(cos(z), sin(z), type="l", 
       main=expression(bold("Circle of Correlations on  ") * bold(list(t[1],t[2]))), 
       ylim=c(-1.1,1.1), xlim=c(-1.2,1.2),
       xlab=expression("PLS-component  " * t[1]), 
       ylab=expression("PLS-component  " * t[2]), 
       cex.main=1, cex.axis=.8, col="grey")
  # adding lines
  abline(h=seq(-1,1,.25), v=seq(-1,1,.25), col="grey", lty=3)
  abline(h=0, v=0, col="grey", lwd=2)
  # variables
  points(cor_tx[,1], cor_tx[,2], pch=20, col=rep("blue",nrow(cor_tx)))
  text(cor_tx[,1], cor_tx[,2], labels=rownames(cor_tx), pos=2, 
       col=rep("blue",nrow(cor_tx)), cex=.8)
  # groups
  points(cor_ty[,1], cor_ty[,2], pch=17, cex=.8, col=rep("red",nrow(cor_ty)))
  text(cor_ty[,1], cor_ty[,2], labels=rownames(cor_ty), pos=2, 
       col=rep("red",nrow(cor_ty)), cex=.8)
  #
  invisible(x)
}
