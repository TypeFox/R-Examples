plot.fmdsd <-
function(x, nscore=1:3, sub.title=NULL, fontsize.points = 1.5, ...)
{
inertia=x$inertia$inertia
coor=x$scores[, -1]
if (max(nscore)>ncol(coor))
  stop("The components of nscore must be smaller than the number of score columns in the x$scores data frame")
group=x$scores[, 1]
ind=combn(nscore, 2)
for (j in 1:ncol(ind)) 
  {i1=ind[1, j]; i2=ind[2, j]
  if (.Device %in% c("null device", "X11", "windows", "quartz"))
    {dev.new()
    }
  par(ps=12);
  plot(coor[,i1],coor[,i2],type="n",
	main="Functional PCA of probability densities",sub=sub.title,
	xlab = paste("PC", i1, " (", inertia[i1], "%)"),ylab=paste("PC", i2, 
         " (", inertia[i2], "%)"), ...)
  par(ps=10);
  text(coor[,i1], coor[,i2], as.character(group), cex=1.5, font=2);
  }
return(invisible(NULL))
}
