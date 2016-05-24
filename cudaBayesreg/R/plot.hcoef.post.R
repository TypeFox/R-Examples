# plot arrays of draws of coefs in hier models
#	3 dimensional arrays:		unit x var x draw
# Original code by P. Rossi 2/07
#
# S3 method to plot arrays of draws of coefs in hier models
#   3 dimensional arrays:  unit x var x draw
# P. Rossi 2/07
#	Modified for cudaBayesreg by A. Ferreira da Silva
#
plot.hcoef.post=function (x, spmname, spm, burnin = trunc(0.1 * R), nsamp=30,   ...) 
{
    X = x
    if (mode(X) == "list") 
        stop("list entered \n Possible Fixup: extract from list \n")
    if (mode(X) != "numeric") 
        stop("Requires numeric argument \n")
    d = dim(X)
    if (length(d) != 3) 
        stop("Requires 3-dim array \n")
    op = par(no.readonly = TRUE)
    on.exit(par(op))
    nunits = d[1]
    nvar = d[2]
    R = d[3]
    if (R < 100) {
        cat("fewer than 100 draws submitted \n")
        return(invisible())
    }
    #
    #	plot posterior distributions of nvar coef for 30 rand units
    #
    len <- length(spm)
    nx <- min(nsamp, len)
    rsam <- sort(sample(spm, nx))
    print(rsam)
    nc <- ceiling(sqrt(nvar))
    nr <- ceiling(nvar/nc)
    par(mfrow = c(nr, nc))
    par(las = 3) # horizontal labeling 
    for (var in 1:nvar) {
        ext = X[rsam, var, (burnin + 1):R]
        ext = data.frame(t(ext))
        colnames(ext) = as.character(rsam)
        out = boxplot(ext, plot = FALSE, ...)
        # !!! change default boxplot.stats based on 1st and 3rd quartiles 
        out$stats = apply(ext, 2, quantile, probs = c(0, 0.05, 
            0.5, 0.95, 1))
        #	bxp(out,xlab="Cross-sectional Unit",title=paste("Coefficients on Var ",
        #   var,sep=""),boxfill="magenta", main=paste(spmname,"var",var),...) 
        bxp(out, xlab = "voxel", title = paste("Coefficients on Var ", 
            var, sep = ""), boxfill = "lightgray", main = paste(spmname, 
            "var", var), boxwex = 0.3, ...)
        if (var == 1) 
            par(ask = dev.interactive())
    }
    par(mfrow = c(1, 1))
    invisible()
}
