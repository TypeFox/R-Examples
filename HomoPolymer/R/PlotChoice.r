PlotChoice<-function(out,otab){
	plotgraph<-function(out,vx,vy,vz=NULL){
		par(mar=c(5,4,4,5)+.1)
		xlim<-c(min(out[,vx]),max(out[,vx]))
		ylim<-c(min(out[,vy]),max(out[,vy]))
		plot(out[,vx],out[,vy[1]],xlim=xlim,ylim=ylim,type='n',
			 xlab=otab$text[ix],ylab=otab$text[iy],col.lab='red')
		grid()
		for(i in 1:length(vy))points(out[,vx],out[,vy[i]],col='red')
		for(i in 1:length(vy))lines(out[,vx],out[,vy[i]],col='red')
		if(!is.null(vz)){
			par(new=TRUE)
			zlim<-c(min(out[,vz]),max(out[,vz]))
			plot(out[,vx],out[,vz[1]],type='n',xaxt="n",yaxt="n",xlab=otab$text[ix],
			ylab=otab$text[iy],,col.lab='red',xlim=xlim,ylim=zlim)
			for(i in 1:length(vz))points(out[,vx],out[,vz[i]],col='blue')
			for(i in 1:length(vz))lines(out[,vx],out[,vz[i]],col='blue')
			axis(4)
			zlab<-otab$text[iz]
			mtext(zlab,side=4,line=3,col='blue')
		}
	}
	var.list<-otab$text
	ans<-inpboxck(c('Variable on x-axis','Secondary y-axis'),var.list,c('-1','FALSE'))
	ix<-as.numeric(ans[[1]])
	vx<-which(names(out)==otab$symbol[ix])
	vz<-ans[[2]]
	ans<-inpboxlist(var.list,'Variables on y main axis')
	iy<-as.numeric(ans[[1]])
	vy<-which(names(out)==otab$symbol[iy])
	if(as.logical(vz)){
	   ans<-inpboxlist(var.list,'Variables on y secondary axis')
	   iz<-as.numeric(ans[[1]])
	   vz<-which(names(out)==otab$symbol[iz])
	   plotgraph(out,vx,vy,vz)
	}else{
	   plotgraph(out,vx,vy)
	}
}

