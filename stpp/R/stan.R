
.listmerge = function (x, y, ...) 
{
### taken from RCurl
    if (length(x) == 0) 
        return(y)
    if (length(y) == 0) 
        return(x)
    i = match(names(y), names(x))
    i = is.na(i)
    if (any(i)) 
        x[names(y)[which(i)]] = y[which(i)]
    x
}


.stan3d.redraw <- function(o) {
  ## switch off redraws
  par3d(skipRedraw=TRUE)
  np=dim(o$xyt)[1]

  ## compute new states
  tin = rep(1,np)
  tin[o$xyt[,3]>(o$t-o$width)]=2
  tin[o$xyt[,3]>o$t]=3

  ## which points have changed state since last time?
  changed = tin != o$xyt[,5]

  ## remove any that changed
  if(any(changed)){
    rgl.pop(id=o$xyt[changed,4])
  }

  ## now add them back in their correct state:
  for(i in (1:np)[changed]){

### should be as simple as this:
###    material3d(o$states[[tin[i]]])
### but setting alpha is causing problems. Bug reported. Hence:
    
    sphereList = list(x=o$xyt[i,1],y=o$xyt[i,2],z=o$xyt[i,3],radius=o$states[[tin[i]]]$radius)
    materialList = o$states[[tin[i]]]
    pList = .listmerge(sphereList,materialList)
    o$xyt[i,4]=do.call(spheres3d,pList)

    o$xyt[i,5]=tin[i]
  }
  ## start drawing again:
  par3d(skipRedraw=FALSE)

  ## update the slider positions
  .store(o)
  
  ## and return the modified object:
  return(o)
}

.rp.stan3d <- function(xyt,tlim,twid,states) {
  t=tlim[1];width=twid
  e=new.env()
  stan.panel  <- rp.control(title="space-time animation",  
                            xyt=xyt, time=tlim[1], width=twid,
                            states=states,
                            e=e
                            )
  Sys.sleep(0.5)
  rp.slider(stan.panel, time, title = "time", from=tlim[1], to=tlim[2], action = .stan3d.redraw,showvalue=TRUE)
  rp.slider(stan.panel, width, title = "window", from=0, to=diff(tlim), action = .stan3d.redraw,showvalue=TRUE)
  rp.button(stan.panel,action=function(p){par3d(FOV=0,userMatrix = rotationMatrix(0, 1,0,0));return(p)},title="reset axes")
  rp.button(stan.panel,action=.store,title="quit",quitbutton=TRUE)
  rp.do(stan.panel, .stan3d.redraw)
  rp.block(stan.panel)
  return(e)
}

.store = function(panel){
 assign("t",panel$time,envir=panel$e)
 assign("width",panel$width,envir=panel$e)
 return(panel)
}

stan <- function(xyt,tlim=range(xyt[,3],na.rm=TRUE),twid=diff(tlim)/20,persist=FALSE,states,bgpoly,bgframe=TRUE,bgimage,bgcol=gray(seq(0,1,len=12)),axes=TRUE){
  require(rgl)
  require(rpanel)

X=(xyt[,1]-min(xyt[,1]))/(diff(range(xyt[,1])))
Y=(xyt[,2]-min(xyt[,2]))/(diff(range(xyt[,2])))
Z=(xyt[,3]-min(xyt[,3]))/(diff(range(xyt[,3])))

tlim=(tlim-min(xyt[,3]))/(diff(range(xyt[,3])))

if(!missing(bgpoly)) bgpoly=cbind((bgpoly[,1]-min(xyt[,1]))/(diff(range(xyt[,1]))),(bgpoly[,2]-min(xyt[,2]))/(diff(range(xyt[,2]))))


XYT=as.3dpoints(X,Y,Z)

  if(missing(states)){
    ## default colouring scheme:
    states=list(
      past=list(col="blue",radius=1/100,alpha=0.5,lit=FALSE),
      present=list(col="red",radius=1/60,alpha=0.5,lit=FALSE),
      ## still-to-come points are invisible (alpha=0)
      future=list(col="yellow",alpha=0.0,radius=1/80,lit=FALSE)
      )
    if(persist){
      states$past=states$present
    }
  }

  maxRadius = max(states$past$radius,states$present$radius,states$future$radius)
  xr = range(XYT[,1],na.rm=TRUE)
  yr = range(XYT[,2],na.rm=TRUE)
  tr = range(XYT[,3],na.rm=TRUE)
  diag=sqrt(diff(xr)^2+diff(yr)^2+diff(tr)^2)

  states$past$radius=states$past$radius*diag
  states$present$radius=states$present$radius*diag
  states$future$radius=states$future$radius*diag
    
  .setPlot(xr[1],xr[2],yr[1],yr[2],tr[1],tr[2],maxRadius,xyt=xyt,axes=axes)

  if(!missing(bgpoly)){
    poly=rbind(bgpoly,bgpoly[1,])
    poly=cbind(poly,min(tr))
    lines3d(poly,size=2.0)
    if(bgframe){
      ci=chull(bgpoly)
      nci=length(ci)
      cpoints=bgpoly[ci,]
      cpoints2 = cpoints[rep(1:nci,rep(2,nci)),]
      cpoints2 = cbind(cpoints2,c(min(tr),max(tr)))
      segments3d(cpoints2)
      poly=cbind(poly[,1:2],max(tr))
      lines3d(poly,size=2.0)
    }
  }

  if(!missing(bgimage)){
    .setBG(bgimage,min(tr),col=bgcol)
  }
  
  XYT=data.frame(XYT[,1:3])
  XYT$id=NA
  ## initially all points will need redrawing:
  XYT$state=-1

  
  ## these points will get redrawn immediately... probably a better way to do this:
  cat("Setting up...")
  npts = dim(XYT)[1]
  if(npts>=100){
    tenths = as.integer(seq(1,npts,len=10))
  }
    
  for(i in 1:(dim(XYT)[1])){
    if(npts>=100){
      tn = (1:10)[tenths == i]
      if(length(tn)>0){
        cat(paste(11-tn,"...",sep=""))
      }
    }
    XYT[i,4]=points3d(XYT[i,1,drop=FALSE],XYT[i,2,drop=FALSE],XYT[i,3,drop=FALSE],alpha=0.0)
  }
  cat("...done\n")

  env = .rp.stan3d(XYT,tlim,twid,states)
  ret=list()
  for(n in ls(env)){
    ret[[n]]=get(n,envir=env)
  }
  return(ret)
}

.ranger <- function(x,margin=0.2){
  lim=range(x,na.rm=TRUE)
  lim=lim+c(-margin,margin)*diff(lim)
  return(lim)
}

.setPlot=function(xmin,xmax,ymin,ymax,tmin,tmax,radius=1/20,xyt,axes){
  require(rgl)
  diag=sqrt((xmax-xmin)^2+(ymax-ymin)^2+(tmax-tmin)^2)
  xr=c(xmax,xmin)+c(radius*(xmax-xmin),-radius*(xmax-xmin))*2
  yr=c(ymax,ymin)+c(radius*(ymax-ymin),-radius*(ymax-ymin))*2
  tr=c(tmax,tmin)+c(radius*(tmax-tmin),-radius*(tmax-tmin))*2
  
  plot3d(xr,yr,tr,type="n",col="red",box=TRUE,axes=FALSE,xlab="x",ylab="y",zlab="t")

	if(axes==TRUE)
	{
	xr2=min(xyt[,1])+range(xr)*(diff(range(xyt[,1])))
 	yr2=min(xyt[,2])+range(yr)*(diff(range(xyt[,2])))
	tr2=min(xyt[,3])+range(tr)*(diff(range(xyt[,3])))
	xl=pretty(xr2,n=5)
	yl=pretty(yr2,n=5)
	tl=pretty(tr2,n=5)
	

	axis3d('x-',at=seq(0,1,length=5),labels=xl,nticks=5)
  	axis3d('y-',at=seq(0,1,length=5),labels=yl,nticks=5)
  	axis3d('z-',at=seq(0,1,length=5),labels=tl,nticks=5)
  	par3d(FOV=0)
  	AR=(xmax-xmin)/(ymax-ymin)
  	aspect3d(AR,1,1)
  	par3d(userMatrix = rotationMatrix(0, 1,0,0))
	}
	else
	{
	axis3d('x-',labels="")
  	axis3d('y-',labels="")
  	axis3d('z-',labels="")
  	par3d(FOV=0)
  	AR=(xmax-xmin)/(ymax-ymin)
  	aspect3d(AR,1,1)
  	par3d(userMatrix = rotationMatrix(0, 1,0,0))
	}
}

.setBG=function(xyz,zplane,col=heat.colors(12)){
  cols = col[as.integer(cut(xyz$z,length(col)))]
  surface3d(xyz$x,xyz$y,rep(zplane,prod(dim(xyz$z))),col=cols,lit=FALSE)
}
