#Plots the polygons in pal imported by get.paldata according
#to the arcs in arc (value returned by get.arcdata)
plotpal<-function(arc,pal, new=TRUE, index=NULL,...)
{

	if(is.null(index))
		index<-1:length(pal[[2]])
	#We only need the lists of arcs
	arc<-arc[[2]]
	pal<-pal[[2]]

	larc<-length(arc)
	arcs<-vector(mode="logical", length=larc)


	for(p in index)
	{
		l<-length(pal[[p]][[1]])
		for(a in 1:l)
		{
			arcs[abs(pal[[p]][[1]][a])]<-TRUE	
		}
	}

	plotarc(list(c(0),arc[arcs]), new=new, ...)

}
