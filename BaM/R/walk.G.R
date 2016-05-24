# Description: 	Plot.walk.g code used to produce figure 9.2
# Usage:	plot.walk.G(walk.mat,sim.rm,X=1,Y=2)



walk <- function(walk.mat,sim.rm,X=1,Y=2)
UseMethod("plot")
walk.G <- function(walk.mat,sim.rm,X=1,Y=2)
{
    plot(walk.mat[1,X],walk.mat[1,Y],type="n",
        xlim=range(walk.mat[,X]),
        ylim=range(walk.mat[,Y]),
	xlab="",ylab="")
    for(i in 1:(nrow(walk.mat)-1))  {
        segments(walk.mat[i,X],walk.mat[i,Y],
                 walk.mat[(i+1),X],walk.mat[i,Y])
        segments(walk.mat[(i+1),X],walk.mat[i,Y],
                 walk.mat[(i+1),X],walk.mat[(i+1),Y])
}
}
