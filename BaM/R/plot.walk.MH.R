# Description: 	Plot.walk.MH code used to produce figure 9.2
# Usage:	plot.walk.MH(walk.mat)


mhplot <- function(walk.mat)
UseMethod("mhplot")
mhplot.walk.MH <- function(walk.mat)  {
    mhplot(walk.mat[1,1],walk.mat[1,2],type="n",
        xlim=round(range(walk.mat[,1])*1.2),
        ylim=round(range(walk.mat[,2])*1.2),
	xlab="",ylab="")
    for(i in 1:(nrow(walk.mat)-1))  {
        segments(walk.mat[i,1],walk.mat[i,2],
                 walk.mat[(i+1),1],walk.mat[(i+1),2])
    }
}	 


