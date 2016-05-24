plot_per_pars <-
function(lca.obj){
  oldpar <- par(no.readonly=T)
  par(mfrow=c(1,2), mar=c(5, 5, 4, 1) + 0.1)
  plot(lca.obj$bx, xlab='Age', ylab=expression(beta[x]^(1)), type='l')
  plot(lca.obj$kt, xlab='Year of birth', ylab=substitute(kappa[t][' '] (a), list(a=lca.obj$adj)))
  title(paste(lca.obj$lab, names(lca.obj)[4], sep=' : '), outer=T, line=-1.5, cex.main=1.5)
  title(lca.obj$model, outer=T, line=-3, cex.sub=0.7, font.sub=3)
  invisible(par(oldpar))
}
