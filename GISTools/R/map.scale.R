`map.scale` <- function (xc,yc,len,units,ndivs,subdiv=1,tcol='black',scol='black',sfcol='black') {
	frame = par("usr")
  l <- len
	tic = (frame[4] - frame[3])/100
	ul = l/ndivs
	for (i in seq(0,ndivs-1,by=2)) rect(xc-l/2+i*ul,yc,xc-l/2+(i+1)*ul,yc+tic/2,border=NA,col=sfcol)
	lines(c(xc-l/2,xc-l/2,xc+l/2,xc+l/2),c(yc+tic,yc,yc,yc+tic),col=scol) 
	lines(c(xc-l/2,xc+l/2),c(yc+tic/2,yc+tic/2),col=scol)
	for (i in 0:ndivs) text(xc-l/2+ul*i,yc-strheight(i*subdiv)*0.7,i*subdiv,col=tcol)
	text(xc,yc-2*strheight(units),units,col=tcol) }

north.arrow <- function (xb,yb,len,lab='NORTH',cex.lab=1,tcol='black',...) {
  s <- len
	arrow.x = c(-1,1,1,1.5,0,-1.5,-1,-1)
	arrow.y = c(0,0,2,2,4,2,2,0) 
	polygon(xb+arrow.x*s,yb+arrow.y*s,...)
	text(xb,yb-strheight(lab,cex=cex.lab),lab,cex=cex.lab,col=tcol)}

