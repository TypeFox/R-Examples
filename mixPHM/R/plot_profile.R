`plot_profile` <-
function (x, method = "mean", type = "b", pch = 19, lty = 1, lwd = 1, col = NA, 
          xlab = "Variables", leglab = NA, ylab = NA, main = NA, legpos = "topright",...)
{
# profile-plot: x-Achse -> Seiten; y-Achse -> Median/Mean der Verweildauern
# method: "mean" or "median"


K <- dim(x$scale)[1]
p <- dim(x$scale)[2]

if (any(is.na(col))) col <- 1:K

if (substitute(method) == "mean") {
  if (is.na(ylab)) ylab <- "Mean of Survival Time"
  if (is.na(main)) main <- "Cluster-Profile \n Mean Survival Times"
  mat <- x$clmean
} else {
  if (is.na(ylab)) ylab <- "Median of Survival Time"
  if (is.na(main)) main <- "Cluster-Profile \n Median Survial Times"
  mat <- x$clmed
}

axlab <- colnames(x$shape)
if (any(is.na(leglab))) leglab <- "Cluster"
if (length(leglab)==1) leglab <- paste(c(leglab),1:K)

#get(getOption("device"))()
matplot(1:p,t(mat),type=type,pch=pch,lty=lty,lwd=lwd,col=col, xaxt="n", xlab=xlab,ylab=ylab,main=main,...)
axis(1,at=1:p,labels=axlab)
legend(legpos,leglab,col=col,lty=lty,lwd=lwd,...)
}

