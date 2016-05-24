star.nominal <-
function (formula, data, xij = NULL, conf.int = FALSE, symmetric = TRUE, 
    pred.coding = "reference", printpvalues = TRUE, test.rel = TRUE, refLevel = 1, 
    maxit = 100, scale = TRUE, nlines = NULL, select = NULL, catstar = TRUE, 
    dist.x = 1, dist.y = 1, dist.cov = 1, dist.cat = 1, xpd = TRUE, main = "", 
    lwd.stars = 1, col.fill = "gray90", col.circle = "black", lwd.circle = 1, 
    lty.circle = "longdash", lty.conf = "dotted", cex.labels = 1, cex.cat = 0.8, 
    xlim = NULL, ylim = NULL)
{ 
if(is.null(xij)){
star.nominalglob(formula, data, conf.int, symmetric, pred.coding, printpvalues, test.rel, 
    refLevel, maxit, scale, nlines, select, dist.x, dist.y, dist.cov, dist.cat, 
    xpd, main, lwd.stars, col.fill, col.circle, lwd.circle, lty.circle, lty.conf, 
    cex.labels, cex.cat, xlim, ylim)
}else{
  if(symmetric){
    warning("Symmetric side constraints are not possible if category-specific covariates are included. 
            The argument 'symmetric' is set 'symmetric = FALSE'.")
    symmetric <- FALSE
  }
star.nominalcat(formula, data, xij, conf.int, symmetric, pred.coding, printpvalues, test.rel, 
    refLevel, maxit, scale, nlines, select, catstar, dist.x, dist.y, dist.cov, 
    dist.cat, xpd, main, lwd.stars, col.fill, col.circle, lwd.circle, lty.circle, 
    lty.conf, cex.labels, cex.cat, xlim, ylim)
}
}
