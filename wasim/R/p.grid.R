p.grid <- function(grid, ...){

myhead <- grid$head
tab <- grid$tab
x=seq(from=myhead[3], by=myhead[5], length.out=myhead[1])
y=seq(from=myhead[4], by=myhead[5], length.out=myhead[2])
filled.contour(x=x,y=y,as.matrix(tab),
	asp=1,
        #color.palette=topo.colors,
	#zlim=c(-20,20),
	plot.axes={axis(1);axis(2)},
	... )
}
