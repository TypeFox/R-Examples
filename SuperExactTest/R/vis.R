#Compute overlap test and visualize intersections between multiple sets
#Author: Minghui Wang
#minghui.wang@mssm.edu
#
plot.msets=function(x,Layout=c('circular','landscape'),degree=NULL,keep.empty.intersections=TRUE,
	sort.by=c('set','size','degree','p-value'),ylim=NULL,
	log.scale=FALSE,x.pos=c(0.05,0.95),y.pos=c(0.025,0.975),yfrac=0.8,color.scale.pos=c(0.85, 0.9),legend.pos=c(0.85,0.25),legend.col=2,legend.text.cex=1,
	color.scale.cex=1,color.scale.title=expression(paste(-Log[10],'(',italic(P),')')),color.on='#2EFE64',color.off='#EEEEEE',
	show.overlap.size=TRUE, show.set.size=TRUE,	track.area.range=0.3,bar.area.range=0.2,origin=if(sort.by[1]=='size'){c(0.45,0.5)}else{c(0.5,0.5)},...){
#keep.empty.intersections, whether to retain empty intersections in the plot
	Layout <- match.arg(Layout)
	if(is.character(color.scale.pos)){
		color.scale.pos=color.scale.pos[1]
	}else if(is.numeric(color.scale.pos)){
		if(length(color.scale.pos)!=2) stop('Invalid color.scale.pos\n')
	}else{
		stop('Invalid color.scale.pos\n')
	}
	sort.by = match.arg(sort.by)
	if(Layout=='circular'){
		return(plot.msets.circular(x=x,degree=degree,keep.empty.intersections=keep.empty.intersections,
		sort.by=sort.by,ylim=ylim,log.scale=log.scale,color.scale.pos=color.scale.pos,legend.pos=legend.pos,
		legend.col=legend.col,legend.text.cex=legend.text.cex,color.scale.cex=color.scale.cex,
		color.scale.title=color.scale.title,color.on=color.on,color.off=color.off,show.overlap.size=show.overlap.size,
		track.area.range=track.area.range,bar.area.range=bar.area.range,origin=origin,...))
	}else if(Layout=='landscape'){
		return(plot.msets.landscape(x=x,degree=degree,keep.empty.intersections=keep.empty.intersections,
		sort.by=sort.by,ylim=ylim,log.scale=log.scale,x.pos=x.pos,y.pos=y.pos,yfrac=yfrac,color.scale.pos=color.scale.pos,color.scale.cex=color.scale.cex,
		color.scale.title=color.scale.title,color.on=color.on,color.off=color.off,show.set.size=show.set.size,...))
	}else{
		stop('Invalid Layout\n')
	}
}
plot.msets.landscape=function(x,degree=NULL,keep.empty.intersections=TRUE,sort.by=c('set','size','degree','p-value'),ylim=NULL,
	log.scale=FALSE,x.pos=c(0.05,1),y.pos=c(0,1),yfrac=0.8,color.scale.pos=c(0.85, 0.9),color.scale.cex=1,color.scale.title=expression(paste(-Log[10],'(',italic(P),')')),
	color.on='#2EFE64',color.off='#EEEEEE',show.set.size=TRUE,...){
	#
	sort.by = match.arg(sort.by)
	Args=list(...)
	cex=ifelse(is.null(Args$cex),par('cex'),Args$cex)
	cex.lab=ifelse(is.null(Args$cex.lab),par('cex.lab'),Args$cex.lab)
	#set color scheme
	if(is.null(Args$heatmapColor)){
		heatmapColor = rev(heat.colors(50)[1:50])
	}else{
		heatmapColor = Args$heatmapColor
	}
	nColors=length(heatmapColor)-1
	params=getPlotParams(x,nColors,degree=degree,keep.empty.intersections=keep.empty.intersections,sort.by=sort.by,ylim=ylim,log.scale=log.scale,Layout='landscape')
	ylabel=params$ylabel
	ylabel0=params$ylabel0
	ylim=params$ylim
	otab=params$otab
	otab[otab>ylim[2]]=ylim[2]
	otab[otab<ylim[1]]=ylim[1]
	otab=otab-ylim[1]
	otab0=params$otab0
	nO=length(otab)
	nSet=length(x$set.sizes) #number of sets
	cid=params$cid
	mlogp=params$mlogp

	#start plotting
	#start a canvas
	grid.newpage()
	vp <- viewport(x=0.5,y=y.pos[1]+(y.pos[2]-y.pos[1])/2, width=1, height=y.pos[2]-y.pos[1])
	pushViewport(vp)
	if(nSet==1){
		grid.circle(x=0.5, y=0.5, r=0.3)
		grid.text(otab0[1],0.5,0.5,gp=gpar(cex=cex*2.5))
		upViewport()
		return(invisible())
	}
	if(nSet==2){
	#to be modified
		grid.circle(x=0.4, y=0.5, r=0.2)
		grid.circle(x=0.6, y=0.5, r=0.2)
		grid.text(otab0['10'],0.3,0.5,gp=gpar(cex=cex*2.5))
		grid.text(otab0['11'],0.5,0.5,gp=gpar(cex=cex*2.5))
		grid.text(otab0['01'],0.7,0.5,gp=gpar(cex=cex*2.5))
		upViewport()
		return(invisible())
	}
	#sub canvas 1
	vp1 <- viewport(x=x.pos[1]+(x.pos[2]-x.pos[1])/2, y=0.6, width=(x.pos[2]-x.pos[1]), height=yfrac)
	pushViewport(vp1)
	
	yLen=1 #height of y axis
	w=1/nO
	h=yLen/(ylim[2]-ylim[1])
	#plot intersections
	for(i in 1:nO){
		posx=w*(i-1)+w/2
	#	grid.text(otab0[i],posx,0.03,rot=45,gp=gpar(cex=cex),just=c('right','top')) #size
		grid.rect(x=posx,y=0.0,width=w*0.8,height=h*otab[i],just=c('center','bottom'),gp=gpar(fill=heatmapColor[cid[i]]))
	}
	upViewport()

	#y axis
	vp1 <- viewport(x=0.05+x.pos[1], y=1-yfrac/2, width=0.1, height=yfrac)
	pushViewport(vp1)
	atVal=(ylabel-ylim[1])*h
	grid.yaxis(at=atVal,label=ylabel0,gp=gpar(cex=cex))
	upViewport()
#	vp1 <- viewport(x=0.00, y=1-yfrac/2, width=0.1, height=yfrac)
#	pushViewport(vp1)
#	grid.text(ylab,1, 0.5, just=c('center'), rot=90, gp=gpar(cex=cex.lab))
#	upViewport()
	
	#color scale
	if((!is.null(x$n)) & (! is.null(mlogp))){
		if(is.character(color.scale.pos)){
			if(color.scale.pos == 'topright'){
				y.vp=0.95
				x.vp=0.9
			}else if (color.scale.pos == 'topleft'){
				y.vp=0.95
				x.vp=0.1
			}else{
				stop('Invalid color.scale.pos\n')
			}
		}else if(is.numeric(color.scale.pos)){
			x.vp=color.scale.pos[1]
			y.vp=color.scale.pos[2]
		}else{
			stop('Invalid color.scale.pos\n')
		}
		vp11 <- viewport(x=x.vp, y=y.vp, width=0.2, height=0.1)
		pushViewport(vp11)
		grid.text(color.scale.title,0.5,0.75,just=c('center','bottom'),gp=gpar(cex=color.scale.cex))
		wc=1/length(heatmapColor)
		for(i in 2:length(heatmapColor)){
			grid.polygon(c((i-1)*wc,i*wc,i*wc,(i-1)*wc),c(0.45,0.45,0.7,0.7),gp=gpar(fill=heatmapColor[i],col=NA))
		}
		t1=0 #floor(min(mlogp,na.rm=T));
		t2=ceiling(max(mlogp,na.rm=T))
		grid.text(t1,0,0.3,just=c('center','top'),gp=gpar(cex=color.scale.cex))
		grid.lines(x = c(wc, wc),y = c(0.45, 0.35))
		grid.text(t2,1,0.3,just=c('center','top'),gp=gpar(cex=color.scale.cex))
		grid.lines(x = c(1-wc, 1-wc),y = c(0.45, 0.35))
		t3=(t1+t2)/2
		if(t3-t1>2){
			grid.text(as.integer(t3),0.5,0.3,just=c('center','top'),gp=gpar(cex=color.scale.cex))
			grid.lines(x = c(0.5, 0.5),y = c(0.45, 0.35))
		}
		upViewport()
	}

	#sub canvas 2, plot matrix
	vp2 <- viewport(x=x.pos[1]+(x.pos[2]-x.pos[1])/2, y=(1-yfrac)/2, width=(x.pos[2]-x.pos[1]), height=1-yfrac)
	pushViewport(vp2)
	h=0.9/nSet
	for(i in 1:nO){
		posx=w*(i-1)+w/2
		for(j in 1:nSet){
			vpJ <- viewport(x=posx, y=0.01+(j-0.5)*h, width=w*0.8, height=h*0.75)
			pushViewport(vpJ)
			if(substr(names(otab[i]),j,j)=='1'){
				grid.circle(0.5,0.5,0.5,gp=gpar(fill=color.on))
			}else{
				grid.circle(0.5,0.5,0.5,gp=gpar(fill=color.off))
			}
			upViewport()
		}
	}
	for(j in 1:nSet){
		grid.text(x$set.names[j], 0.1*w, 0.015+(j-0.5)*h,just=c('right','center'),gp=gpar(cex=cex))
		if(show.set.size==TRUE) grid.text(x$set.sizes[j], 1, 0.015+(j-0.5)*h,just=c('left','center'),gp=gpar(cex=cex))
	}
	upViewport()
	upViewport()
	return(invisible())
}
plot.msets.circular=function(x,degree=NULL,keep.empty.intersections=TRUE,sort.by=c('set','size','degree','p-value'),
	ylim=NULL,log.scale=FALSE,color.scale.pos=c(0.85, 0.9), legend.pos=c(0.85,0.25),legend.col=2,legend.text.cex=1,
	color.scale.cex=1,color.scale.title=expression(paste(-Log[10],'(',italic(P),')')),color.on='#2EFE64',color.off='#EEEEEE',
	show.overlap.size=TRUE,track.area.range=0.3,bar.area.range=0.2,origin=if(sort.by[1]=='size'){c(0.45,0.5)}else{c(0.5,0.5)},...){
	#
	if(is.character(legend.pos)){
		legend.pos <- legend.pos[1]
	}else if(is.numeric(legend.pos)){
		if(length(legend.pos)!=2) stop('Invalid legend.pos\n')
	}else{
		stop('Invalid legend.pos\n')
	}
	sort.by = match.arg(sort.by)
	Args=list(...)
	cex=ifelse(is.null(Args$cex),0.8,Args$cex)
	#set color scheme
	if(is.null(Args$heatmapColor)){
		heatmapColor = rev(heat.colors(50)[1:50])
	}else{
		heatmapColor = Args$heatmapColor
	}
	nColors=length(heatmapColor)-1
	params=getPlotParams(x,nColors,degree=degree,keep.empty.intersections=keep.empty.intersections,sort.by=sort.by,ylim=ylim,log.scale=log.scale)
	ylim=params$ylim
	otab=params$otab
	otab[otab>ylim[2]]=ylim[2]
	otab[otab<ylim[1]]=ylim[1]
	otab=otab-ylim[1]
	otab0=params$otab0
	nO=length(otab)
	nSet=length(x$set.sizes) #number of sets
	cid=params$cid
	mlogp=params$mlogp
	radial2deg=180/pi

	# set graph layout parameters
	track.area.range=as.numeric(track.area.range)
	bar.area.range=as.numeric(bar.area.range)
	if(track.area.range  < 0) stop('Invalid argument "trak.area.range"\n')
	if(bar.area.range  < 0) stop('Invalid argument "bar.area.range"\n')
	if(track.area.range+bar.area.range>0.5) stop('Sum of track.area.range and bar.area.range should not be larger than 0.5\n')
	width.sets=track.area.range
	width.intersections=bar.area.range
	track.offset=ifelse(is.null(Args$phantom.traks),2,as.integer(Args$phantom.traks)) #number of phantom tracks in the middle
	track.width=width.sets/(nSet+track.offset)
	bar.width.unit=width.intersections/(ylim[2]-ylim[1])
	gap.within.track=ifelse(is.null(Args$gap.within.track),0.1,Args$gap.within.track) #ratio of gap width over block width on the same track
	gap.between.track=ifelse(is.null(Args$gap.between.track),0.1,Args$gap.between.track) #ratio of gap width over track width

	#start a canvas
	grid.newpage()

	#Plot tracks
	
	vp1 <- viewport(x=origin[1],y=origin[2],width=0.95, height=0.95)
	pushViewport(vp1)
	degreeUnit=2*pi/nO
	degreeStart=(c(1:nO)-1)*degreeUnit
	degreeEnd=(c(1:nO))*degreeUnit
	degree.gap=max(degreeUnit*gap.within.track,2*pi/360)

	fill.col=rep(c('#dddddd','#999999'),length.out=nSet)
	#plot tracks
	#for(i in 1:nSet){
	#	XY1=sapply(seq(0,2*pi,length.out=360), function(deg) getXY(origin,(i+track.offset-1)*track.width,deg))
	#	XY2=sapply(seq(0,2*pi,length.out=360), function(deg) getXY(origin,(i+track.offset)*track.width-track.width*gap.between.track,deg))
	#	pos.x=c(XY1[1,],rev(XY2[1,]))
	#	pos.y=c(XY1[2,],rev(XY2[2,]))
	#	grid.polygon(pos.x, pos.y,gp=gpar(col=1))
	#}
	#plot overlap (intersection)
	for(i in 1:nO){
		which.set=strsplit(names(otab)[i],'')[[1]]=='1'
		for(j in 1:nSet){
			XY1=sapply(seq(degreeStart[i],degreeEnd[i]-degree.gap,length.out=40), function(deg) getXY(origin,(j+track.offset-1)*track.width,deg))
			XY2=sapply(seq(degreeStart[i],degreeEnd[i]-degree.gap,length.out=40), function(deg) getXY(origin,(j+track.offset)*track.width-track.width*gap.between.track,deg))
			pos.x <- c(XY1[1,],rev(XY2[1,]))
			pos.y <- c(XY1[2,],rev(XY2[2,]))
			grid.polygon(pos.x, pos.y,gp=gpar(col = 'black',fill=ifelse(which.set[j],color.on,color.off))) #
		}
		#bar plot intersection size
		XY1=sapply(seq(degreeStart[i],degreeEnd[i]-degree.gap,length.out=40), function(deg) getXY(origin,width.sets,deg))
		XY2=sapply(seq(degreeStart[i],degreeEnd[i]-degree.gap,length.out=40), function(deg) getXY(origin,width.sets+bar.width.unit*otab[i],deg))
		pos.x=c(XY1[1,],rev(XY2[1,]))
		pos.y=c(XY1[2,],rev(XY2[2,]))
		grid.polygon(pos.x, pos.y,gp=gpar(fill = heatmapColor[cid[i]],col=1)) #bar height is proportional to the intersection size
		#text intersection size
		if(show.overlap.size==TRUE){
			XY3=sapply(seq(degreeStart[i],degreeEnd[i]-degree.gap,length.out=4), function(deg) getXY(origin,width.sets+0.005+bar.width.unit*otab[i],deg))
			grid.text(otab0[i],mean(XY3[1,]),mean(XY3[2,]),rot=degreeStart[i]*radial2deg,just='left',gp=gpar(cex=cex))
		}
	}
	#track number (numbering the legends)
	for(j in 1:nSet){
		XY1=sapply(seq(degreeStart[1],degreeEnd[1]-degree.gap,length.out=40), function(deg) getXY(origin,(j+track.offset-1)*track.width,deg))
		XY2=sapply(seq(degreeStart[1],degreeEnd[1]-degree.gap,length.out=40), function(deg) getXY(origin,(j+track.offset)*track.width,deg))
		grid.text(j,(XY1[1,1]+XY2[1,1])/2,mean(XY1[2,]),just=c('center'),gp=gpar(cex=cex))
	}
	upViewport()
	#color scale
	if((!is.null(x$n)) & (! is.null(mlogp))){
		if(is.character(color.scale.pos)){
			if(color.scale.pos == 'topright'){
				y.vp=0.90
				x.vp=0.85
			}else if(color.scale.pos == 'bottomright'){
				y.vp=0.15
				x.vp=0.85
			}else if (color.scale.pos == 'topleft'){
				y.vp=0.90
				x.vp=0.15
			}else if(color.scale.pos == 'bottomleft'){
				y.vp=0.15
				x.vp=0.15
			}else{
				stop('Invalid color.scale.pos\n')
			}
		}else if(is.numeric(color.scale.pos)){
			x.vp=color.scale.pos[1]
			y.vp=color.scale.pos[2]
		}else{
			stop('Invalid color.scale.pos\n')
		}
		vp11 <- viewport(x=x.vp, y=y.vp, width=0.2, height=0.1)
		pushViewport(vp11)
		grid.text(color.scale.title,0.5,0.75,just=c('center','bottom'),gp=gpar(cex=color.scale.cex))
		wc=1/length(heatmapColor)
		for(i in 2:length(heatmapColor)){
			grid.polygon(c((i-1)*wc,i*wc,i*wc,(i-1)*wc),c(0.45,0.45,0.7,0.7),gp=gpar(fill=heatmapColor[i],col=NA))
		}
		t1=0 #floor(min(mlogp,na.rm=T));
		t2=ceiling(max(mlogp,na.rm=T))
		grid.text(t1,0,0.3,just=c('center','top'),gp=gpar(cex=color.scale.cex))
		grid.lines(x = c(wc, wc),y = c(0.45, 0.35))
		grid.text(t2,1,0.3,just=c('center','top'),gp=gpar(cex=color.scale.cex))
		grid.lines(x = c(1-wc/2, 1-wc/2),y = c(0.45, 0.35))
		t3=(t1+t2)/2
		if(t3-t1>2){
			grid.text(as.integer(t3),0.5,0.3,just=c('center','top'),gp=gpar(cex=color.scale.cex))
			grid.lines(x = c(0.5, 0.5),y = c(0.45, 0.35))
		}
		upViewport()
	}

	#plot legend
	if(is.character(legend.pos)){
		if(legend.pos == 'topright'){
			y.vp=0.80
			x.vp=0.85
		}else if(legend.pos == 'bottomright'){
			y.vp=0.25
			x.vp=0.85
		}else if (legend.pos == 'topleft'){
			y.vp=0.80
			x.vp=0.15
		}else if(legend.pos == 'bottomleft'){
			y.vp=0.25
			x.vp=0.15
		}else{
			stop('Invalid legend.pos\n')
		}
	}else if(is.numeric(legend.pos)){
		x.vp=legend.pos[1]
		y.vp=legend.pos[2]
	}else{
		stop('Invalid legend.pos\n')
	}
	legend.box.width=0.2
	legend.box.height=0.2
	vp2 <- viewport(x=x.vp, y=y.vp, width=legend.box.width, height=legend.box.height)
	pushViewport(vp2)
	Legend <- frameGrob()
	Legend <- packGrob(Legend, grid.legend(paste(1:nSet,x$set.names,sep=': '),ncol=legend.col,
		vgap=ifelse(is.null(Args$legend.vgap),0.1*legend.text.cex,Args$legend.vgap),hgap=ifelse(is.null(Args$legend.hgap),0.2*legend.text.cex,Args$legend.hgap),
		gp=gpar(cex=legend.text.cex), draw = FALSE), height = unit(1, "null"),side = ifelse(x.vp>0.5,"right",'left'))
	grid.draw(Legend)
	upViewport()
	return(invisible())
}
getXY=function(origin,radius,degree){
	X=radius*cos(degree)
	Y=radius*sin(degree)
	origin+c(X,Y)
}
getPlotParams=function(x,nColors=50,degree=NULL,keep.empty.intersections=TRUE,sort.by=c('set','size','degree','p-value'),ylim=NULL,log.scale=FALSE,Layout=c('circular','landscape')){
	Layout=match.arg(Layout)
	sort.by = match.arg(sort.by)
	otab=x$overlap.sizes
	if(sort.by=='set'){
		otab.order=order(names(otab))
		otab=otab[otab.order]
	}else if(sort.by=='size'){
		otab.order=order(otab,decreasing=TRUE)
		otab=otab[otab.order]
	}else if(sort.by=='degree'){
		otab.order=order(sapply(names(otab),function(d) countCharOccurrences('1',d))) #order(sapply(strsplit(names(otab),''),function(d) sum(d=='1')))
		otab=otab[otab.order]
	}else if(sort.by=='p-value'){
		otab.order=order(x$P.value)
		otab=otab[otab.order]
	}else{
		stop('Invalid sort.by argument\n')
	}
	if(keep.empty.intersections==FALSE){
		otab.kept=otab>0
		otab=otab[otab.kept]
	}else{
		otab.kept=rep(T,length(otab))
	}
	if(!is.null(degree)){
		kpt=sapply(names(otab),function(d) countCharOccurrences('1',d)) %in% degree #sapply(strsplit(names(otab),''),function(d) sum(d == '1') %in% degree)
		if(sum(kpt)<2) stop('Too few items left for plotting after applying degree filter\n')
		otab.kept[otab.kept]=kpt
		otab=otab[kpt]
	}
	otab0=otab
	if(is.null(ylim)){
		ylabel=axisTicks(c(0,max(otab0,na.rm=T)),log=FALSE) #for landscape layout
		ylim=c(0,max(otab,na.rm=T))
		if(Layout=='landscape') ylim[2]=max(ylabel)
	}else{
		ylabel=axisTicks(ylim[1:2],log=FALSE)
	}
	if(ylim[2]<ylim[1]) stop("Invalid ylim\n")
	ylabel0=ylabel
	if(log.scale==TRUE){
		if(any(ylim < 0)) stop("ylim can't be negative when log.scale is TRUE\n")
		otab=log(otab+1);
		ylim=log(ylim+1)
		ylabel=log(ylabel+1)
	}
	nO=length(otab)
	cid=rep(1,nO) #color gradient id
	mlogp=NULL
	if((!is.null(x$n)) & nO>0){
		mlogp=-log(x$P.value,base=10)
		mlogp[is.na(mlogp)]=1.0e-10
		mlogp=mlogp[otab.order]
		mlogp=mlogp[otab.kept]
		mlogp[mlogp == Inf] = max(320,mlogp[mlogp < Inf],na.rm=T)
		cid=ceiling(nColors*mlogp/max(c(mlogp,1e-10),na.rm=T))
		cid[cid<1]=1
		if(all(is.na(mlogp))) mlogp=NULL
	}
	return(list(otab=otab,otab0=otab0,cid=cid,mlogp=mlogp,ylim=ylim,ylabel=ylabel,ylabel0=ylabel0))
}
axisLen=function(Max,n.ticks=5){
	u=ceiling(Max/n.ticks)
	u*n.ticks
}
