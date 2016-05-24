magbar <-
function(position='topright',range=c(0,1),orient='v',log=FALSE,col=hsv(h=seq(2/3,0,len=100)),scale=c(1/4,1/20),inset=1/40,labN=5,title='',titleshift=0,centrealign='rb',clip=''){
usercoord=par()$usr
xlogcheck=FALSE;ylogcheck=FALSE
if(par()$xlog){par(xlog=FALSE);par(usr=c(log10(par()$usr[1:2]),par()$usr[3:4]));xlogcheck=TRUE}
if(par()$ylog){par(ylog=FALSE);par(usr=c(par()$usr[1:2],log10(par()$usr[3:4])));ylogcheck=TRUE}
                 
xlo=usercoord[1];xhi=usercoord[2];ylo=usercoord[3];yhi=usercoord[4]
xdiff=xhi-xlo;ydiff=yhi-ylo

if(orient=='h'){
	xl=xlo+xdiff/2-xdiff*scale[1]/2
	yb=ylo+ydiff/2-ydiff*scale[2]/2
	xr=xlo+xdiff/2+xdiff*scale[1]/2
	yt=ylo+ydiff/2+ydiff*scale[2]/2
	align=centrealign
	if(length(grep('bottom',position))>0){align='lt';yb=ylo+ydiff*inset;yt=ylo+ydiff*inset+ydiff*scale[2]}
	if(length(grep('top',position))>0){align='rb';yb=yhi-ydiff*inset-ydiff*scale[2];yt=yhi-ydiff*inset}
	if(length(grep('left',position))>0){xl=xlo+xdiff*inset;xr=xlo+xdiff*inset+xdiff*scale[1]}
	if(length(grep('right',position))>0){xl=xhi-xdiff*inset-xdiff*scale[1];xr=xhi-xdiff*inset}
}

if(orient=='v'){
	xl=xlo+xdiff/2-xdiff*scale[2]/2
	yb=ylo+ydiff/2-ydiff*scale[1]/2
	xr=xlo+xdiff/2+xdiff*scale[2]/2
	yt=ylo+ydiff/2+ydiff*scale[1]/2
	align=centrealign
	if(length(grep('bottom',position))>0){yb=ylo+ydiff*inset;yt=ylo+ydiff*inset+ydiff*scale[1]}
	if(length(grep('top',position))>0){yb=yhi-ydiff*inset-ydiff*scale[1];yt=yhi-ydiff*inset}
	if(length(grep('left',position))>0){align='rb';xl=xlo+xdiff*inset;xr=xlo+xdiff*inset+xdiff*scale[2]}
	if(length(grep('right',position))>0){align='lt';xl=xhi-xdiff*inset-xdiff*scale[2];xr=xhi-xdiff*inset}
}

rangetemp=range
legend=maglab(rangetemp,labN,log=log,trim=F)
if(log){
  rangetemp=log10(range)
  legend$labat=log10(legend$labat)
}

roughNscale=(max(legend$labat)-min(legend$labat))/(rangetemp[2]-rangetemp[1])
if(clip=='bg'){clip='NA'}

colremap=magmap(data=seq(min(legend$labat),max(legend$labat),length=length(col)*roughNscale),lo=rangetemp[1],hi=rangetemp[2],type='num',range=c(1,length(col)),clip=clip)$map
col=col[round(colremap,digits=0)]

if(orient=='v'){color.legend(xl,yb,xr,yt,legend=legend$exp,rect.col=col,align=align,gradient='y')}
if(orient=='h'){color.legend(xl,yb,xr,yt,legend=legend$exp,rect.col=col,align=align,gradient='x')}

if(orient=='v' & align=='lt'){text(xl-(1+titleshift)*xdiff/20,(yt+yb)/2,labels=title,adj=c(0.5,0.5),srt=90)}
if(orient=='v' & align=='rb'){text(xr+(1+titleshift)*xdiff/20,(yt+yb)/2,labels=title,adj=c(0.5,0.5),srt=-90)}
if(orient=='h' & align=='lt'){text((xl+xr)/2,yt+(1+titleshift)*ydiff/20,labels=title,adj=c(0.5,0.5),srt=0)}
if(orient=='h' & align=='rb'){text((xl+xr)/2,yb-(1+titleshift)*ydiff/20,labels=title,adj=c(0.5,0.5),srt=0)}

par(xlog=xlogcheck)
par(ylog=ylogcheck)
par(usr=usercoord)
}

