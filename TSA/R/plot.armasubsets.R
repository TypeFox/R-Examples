`plot.armasubsets` <-
function (x, labels = obj$xnames, main = NULL, scale = c("BIC","AICc","AIC", 
    "Cp", "adjR2", "R2"), col = gray(c(seq(0.4, .7, length = 10),.9)),draw.grid=TRUE, 
    axis.at.3=TRUE,...) 
{
# date March 2, 2006
# Programmed by Kung-Sik Chan
#
# modified from plot.subsets with two new options
# draw.grid if it is true (default), gray grid lines are superimposed on the figure.
# axis.at.3 if it is true (default), the x-labels are drawn on the upper 
# horizontal axis.
# scale determines which model selection criteria will be used.
#
    obj <- x
    lsum <- summary.armasubsets(obj)
    par(mar = c(2, 5, 7, 3) + 0.1)
    nmodels <- length(lsum$rsq)
    np <- obj$np
    propscale <- FALSE
    sscale <- pmatch(scale[1], c("BIC","AICc","AIC", "Cp", "adjR2", "R2"), 
        nomatch = 0)
    if (sscale == 0) 
        stop(paste("Unrecognised scale=", scale))
    if (propscale) 
        stop(paste("Proportional scaling only for probabilities"))
    yscale <- switch(sscale,lsum$BIC, lsum$AICc, lsum$AIC,  lsum$Cp, lsum$adjR2, lsum$rsq)
    up <- switch(sscale, -1,-1,-1, -1, 1, 1)
    index <- order(yscale * up)
    colorscale <- switch(sscale, yscale, yscale, -log(pmax(yscale, 
        1e-04)), -log(pmax(yscale, 1e-04)))
    image(z = t(ifelse(lsum$which[index, ], colorscale[index], 
        NA + max(colorscale) * 1.5)), xaxt = "n", yaxt = "n", 
        x = (1:np), y = 1:nmodels, xlab = "", ylab = scale[1], 
        col = col)
if (draw.grid) {
for (i in seq(.5,(np-.5))) abline(h=i,col="gray")
for( j in (seq(labels)+.5)) abline(v=j,col="gray")
}
    laspar <- par("las")
    on.exit(par(las = laspar))
    par(las = 2)
for(i in seq(labels)) {
if(charmatch('`',labels[i],nomatch=0)==1) labels[i]=substr(labels[i],2,nchar(labels[i])-1)
}
    if (axis.at.3) axis(3, at = 1:np, labels = labels) else
    axis(1, at = 1:np, labels = labels)
    axis(2, at = 1:nmodels, labels = signif(yscale[index], 2)) 
    if (!is.null(main)) 
        title(main = main)
    box()
    invisible(NULL)
}

