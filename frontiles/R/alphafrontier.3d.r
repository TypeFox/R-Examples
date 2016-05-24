alphafrontier.3d<-function(xobs, yobs, type="output",alpha=0.95, digits=4, box.leg=TRUE,
palette=heat_hcl, rgl=FALSE, n.class=NULL,  ...)
{
# initialisation
m <- match.call(expand.dots = FALSE)
n<-nrow(xobs)
p<-ncol(xobs)

if(n!=nrow(yobs)) return("xobs and yobs have not the same number of observations")
if(p!=2 | ncol(yobs)!=1) return("xobs must have 2 columns and yobs 1 column")
match.arg(type,c("output","input"))

xobs1<-xobs[,1]
xobs2<-xobs[,2]
 
if(type=="input")
{
yobs.sort<-sort(yobs)
n.k<-0  # discretisation foot
xeval<-matrix(0,2*n.k+length(yobs),2)

k<-1
interval<-yobs.sort[k]
 yeval<-rep(yobs.sort[k],length(yobs))

 res<-alphascore(xobs=xobs,yobs=yobs, xeval=xobs, yeval=as.matrix(yeval),alpha=alpha)
 x.alqf.2<-xeval*res$input

 # For each value of yobs, we verify if the frontier is the same with the previous
 # value of yobs
 
for (k in 2:n)
{
 yeval<-rep(yobs.sort[k],length(yobs))

 res<-alphascore(xobs=xobs,yobs=yobs, xeval=xobs, yeval=as.matrix(yeval),alpha=alpha)
 x.alqf<-xeval*res$input

  if(!all(x.alqf==x.alqf.2))
  {x.alqf.2<-x.alqf
   interval<-c(interval,yobs.sort[k])
  }
}

    n.couleur<-length(interval)
    type.col<-rev(palette(n.couleur))

if((box.leg==TRUE)&(n.couleur>1))
  {layout(matrix(c(1,1,2,2),2,2), widths = c(1,4), heights =c(1,1))
   barplot(rep(1,n.couleur),horiz=TRUE,col=type.col,space=0,axes=FALSE)
   axis(2,0:(n.couleur-1),interval,las=1,cex=0.6)
  }

##############################
##############################
n.k<-100  # discretisation foot
xeval<-matrix(0,2*n.k+length(yobs),2)

# we prepare the bound
bound.x1<-c(min(xobs1)-diff(range(xobs1))/20,max(xobs1)+diff(range(xobs1))/20)
bound.x2<-c(min(xobs2)-diff(range(xobs2))/20,max(xobs2)+diff(range(xobs2))/20)

plot(xobs,type='n')#,xlim=c(300,500),ylim=c(4000,5000))
 for (k in 1:n.couleur)
 {
  yeval<-c(rep(interval[k],2*n.k+length(yobs)))
  xeval[,1]<-c(seq(min(xobs[,1]),max(xobs[,1]),length.out=n.k),seq(mean(xobs[,1]),mean(xobs[,1]),length.out=n.k),xobs[,1])
  xeval[,2]<-c(seq(mean(xobs[,2]),mean(xobs[,2]),length.out=n.k),seq(min(xobs[,2]),max(xobs[,2]),length.out=n.k),xobs[,2])

  res<-alphascore(xobs=xobs,yobs=yobs, xeval=xeval, yeval=as.matrix(yeval),alpha=1)
  x.alqf<-xeval*res$input
  x.alqf<-rbind(x.alqf,rbind(c(bound.x1[2],min(x.alqf[,2])),c(min(x.alqf[,1]),bound.x2[2])))

  x.alqf[,1]<-sort(x.alqf[,1])
  x.alqf[,2]<-sort(x.alqf[,2],decreasing=TRUE)
  lines(x.alqf,col=type.col[k])
 }
 }
 else     ### CASE OUTPUT DIRECTION
 {    # we prepare the bound
bound.x1<-c(min(xobs1)-diff(range(xobs1))/20,max(xobs1)+diff(range(xobs1))/20)
bound.x2<-c(min(xobs2)-diff(range(xobs2))/20,max(xobs2)+diff(range(xobs2))/20)

# we order x1 and x2 to have a grid  of size (n)x(n)
x1.order<-sort(xobs1)
x2.order<-sort(xobs2)

# we "vectorize" the grid
xpred1<-rep(x1.order,n)
xpred2<-rep(x2.order,each=n)
xpred=cbind(xpred1,xpred2)

# we add to x1.ord and x2.ord the bound max
x1.order<-c(x1.order,bound.x1[2]+diff(bound.x1)/3)
x2.order<-c(x2.order,bound.x2[2]+diff(bound.x2)/3)

n2<-nrow(xpred)

# 1- Computation of the alpha-quantile score
  yeval<-matrix(mean(yobs),n2)

# Call to the C.code       
res<-.C("alpha3d",as.integer(n), as.integer(n2), as.double(t(xobs)), as.double(t(yobs)),
as.double(t(xpred)), as.double(t(yeval)), lambda=as.double(matrix(0,n2,1)),res1=as.double(matrix(0,n,1)), 
as.double(alpha) ,PACKAGE="frontiles")

  
 alqf.score<-res$lambda
 alqf.score[alqf.score==-1]<-0
    
  y.frontier<-alqf.score*mean(yobs)  # score computed on (xobs,yobs)
  
   
  # unique values
  y.fr.unique<-unique(na.omit(round(y.frontier,digits)))
  n.fdh<-length(y.fr.unique)
  
 # we divide the frontier into n intervals (use of package classInt)
 # with a special color for each class
 
    n.couleur<-ifelse(is.null(n.class),n.fdh, n.class)
    type.col<-rev(palette(n.couleur))
    class.var<-classIntervals(y.fr.unique, n=n.couleur)$brks
  

 # choice of representation 2-d with colors or 3-d
if(!rgl)
 {if((box.leg==TRUE)&(n.couleur>1))
  {layout(matrix(c(1,1,2,2),2,2), widths = c(1,4), heights =c(1,1))
   barplot(rep(1,n.couleur),horiz=TRUE,col=type.col,space=0,axes=FALSE)
   if(n.couleur==n.fdh)
    {axis(2,0.5:n.couleur,sort(y.fr.unique),las=1,cex=0.6)}
    else
    {axis(2,1:(n.couleur+1),class.var,las=1,cex=0.6)}  
   }

 #we define a SpatialPolygonsDataFrame
 pol<-vector("list", n^2)

 for (k in 1:n)
  { ypol=c(x2.order[k],x2.order[k],x2.order[k+1],x2.order[k+1],x2.order[k])
  for (j in 1:n)
   {xpol=c(x1.order[j],x1.order[j+1],x1.order[j+1],x1.order[j],x1.order[j])
    pol[[n*(k-1)+j]] = Polygons(list(Polygon(cbind(xpol,ypol))), ID=paste("x",n*(k-1)+j,sep=""))
   }
  }

 sr=SpatialPolygons(pol)
 srdf=SpatialPolygonsDataFrame(sr, data.frame(y.frontier, row.names=paste("x",1:n^2,sep="")))
 
 plot(coordinates(srdf),type="n",...)
 plot(srdf,col=type.col[findInterval(srdf@data[,1], class.var, all.inside=TRUE)],
  border=type.col[findInterval(srdf@data[,1], class.var, all.inside=TRUE)],add=TRUE)
 }
 else
 {
  z<-matrix(y.frontier,n,n)
  main<-ifelse(!is.na(pmatch("main",names(m$...))), eval(m$...$main),'')
  xlab<-ifelse(!is.na(pmatch("xlab",names(m$...))), eval(m$...$xlab),'')     
  ylab<-ifelse(!is.na(pmatch("ylab",names(m$...))), eval(m$...$ylab),'')
  zlab<-ifelse(!is.na(pmatch("zlab",names(m$...))),eval(m$...$zlab),'')
  sub<-ifelse(!is.na(pmatch("sub",names(m$...))),eval(m$...$sub),'')
   for (i in 1:n)
    {
    for (j in 1:n)
      {
      binplot.3d(c(x1.order[i],x1.order[i+1]),c(x2.order[j],x2.order[j+1]),
                 z[i,j],alpha=1,topcol=type.col[findInterval(z[i,j], class.var, all.inside=TRUE)])
      }
      aspect3d(1,1,1)
    }
      axes3d(c('x','y','z'))
      title3d(main,sub,xlab,zlab,ylab)
 }
 }
}

 ##### Required functions 'binplot' and 'hist3d':
 
binplot.3d<-function(x,y,z,alpha=1,topcol="#ff0000",sidecol="#aaaaaa")
{
  save <- par3d(skipRedraw=TRUE)
  on.exit(par3d(save))
    
  x1<-c(rep(c(x[1],x[2],x[2],x[1]),3),rep(x[1],4),rep(x[2],4))
  z1<-c(rep(0,4),rep(c(0,0,z,z),4))
  y1<-c(y[1],y[1],y[2],y[2],rep(y[1],4),rep(y[2],4),rep(c(y[1],y[2],y[2],y[1]),2))
  x2<-c(rep(c(x[1],x[1],x[2],x[2]),2),rep(c(x[1],x[2],rep(x[1],3),rep(x[2],3)),2))
  z2<-c(rep(c(0,z),4),rep(0,8),rep(z,8) )
  y2<-c(rep(y[1],4),rep(y[2],4),rep(c(rep(y[1],3),rep(y[2],3),y[1],y[2]),2) )
  rgl.quads(x1,z1,y1,col=rep(sidecol,each=4),alpha=alpha)
  rgl.quads(c(x[1],x[2],x[2],x[1]),rep(z,4),c(y[1],y[1],y[2],y[2]),
              col=rep(topcol,each=4),alpha=1) 
  rgl.lines(x2,z2,y2,col="#000000")
}

