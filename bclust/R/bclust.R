bclust<-function(x,rep.id=1:nrow(x),effect.family="gaussian",
var.select=TRUE,transformed.par,labels=NULL)
{
if (missing(x)){stop("data are not specified")}
if (!is.matrix(x)){stop("data must be a matrix")}
if (missing(transformed.par)){stop("transformed.par is missing")}
if (!is.vector(transformed.par)){stop("transformed.par must be a vector")}
if (!(length(rep.id)==nrow(x))){stop("rep.id mismatches")}
if ((sum(is.na(x))>0) | (sum(is.na(rep.id))>0) | 
(sum(is.na(transformed.par))>0)) stop("NA is not allowed")
if(!(effect.family%in%c("gaussian","alaplace"))) {stop(
"in effect.family just the gaussian and the alaplace families are allowed")}
id.order<-order(rep.id)
repno<-as.vector(table(rep.id[id.order]))
if (is.null(labels)){if (is.null(rownames(x))){labels<-paste(1:length(repno))} else  
{if (length(repno)==nrow(x)){labels<-rownames(x)};labels<-labels[id.order]} } else 
{if (!(length(labels)==length(repno))){warning("label length does not match with the data")}}

y<-x[id.order,]
if (is.null(colnames(y))) {colnames(y)<-paste(1:ncol(y))}

bclustG<-function(y,repno,hyperparameters)
{
out<-.C("RfastbclustG",PACKAGE="bclust",y=as.double(y),nrowy=as.integer(nrow(y)),ncoly=as.integer(ncol(y)),repno=as.double(repno),
nrepno=as.integer(length(repno)),theta=as.double(hyperparameters),merge=as.double(rep(0,2*(length(repno)-1))),
height=as.double(rep(0,(length(repno)-1))))
mytree<-list(merge=matrix(out$merge,ncol=2,byrow=TRUE),height=out$height,order=nocrossing.order(matrix(out$merge,ncol=2,byrow=TRUE))
)
return(mytree)
}

bclustvsG<-function(y,repno,hyperparameters)
{
out<-.C("RfastbclustvsG",PACKAGE="bclust",y=as.double(y),nrowy=as.integer(nrow(y)),ncoly=as.integer(ncol(y)),repno=as.double(repno),
nrepno=as.integer(length(repno)),theta=as.double(hyperparameters),merge=as.double(rep(0,2*(length(repno)-1))),
height=as.double(rep(0,(length(repno)-1))))
mytree<-list(merge=matrix(out$merge,ncol=2,byrow=TRUE),height=out$height,order=nocrossing.order(matrix(out$merge,ncol=2,byrow=TRUE))
)
return(mytree)
}

bclustAL<-function(y,repno,hyperparameters)
{
out<-.C("RfastbclustAL",PACKAGE="bclust",y=as.double(y),nrowy=as.integer(nrow(y)),ncoly=as.integer(ncol(y)),repno=as.double(repno),
nrepno=as.integer(length(repno)),theta=as.double(hyperparameters),merge=as.double(rep(0,2*(length(repno)-1))),
height=as.double(rep(0,(length(repno)-1))))
mytree<-list(merge=matrix(out$merge,ncol=2,byrow=TRUE),height=out$height,order=nocrossing.order(matrix(out$merge,ncol=2,byrow=TRUE))
)
return(mytree)
}


bclustvsAL<-function(y,repno,hyperparameters)
{
out<-.C("RfastbclustvsAL",PACKAGE="bclust",y=as.double(y),nrowy=as.integer(nrow(y)),ncoly=as.integer(ncol(y)),repno=as.double(repno),
nrepno=as.integer(length(repno)),theta=as.double(hyperparameters),merge=as.double(rep(0,2*(length(repno)-1))),
height=as.double(rep(0,(length(repno)-1))))
mytree<-list(merge=matrix(out$merge,ncol=2,byrow=TRUE),
height=out$height,order=nocrossing.order(matrix(out$merge,ncol=2,byrow=TRUE)))
return(mytree)
}

nocrossing.order<-function(mergemat)
{
	myorder<-matrix(NA,nrow(mergemat),nrow(mergemat)+1)
			na.add <- function(vec)
	{
		return(c(vec,rep(NA,nrow(mergemat)+1-length(vec))))
	}

	myorder[1,]<-na.add(mergemat[1,])
			i<-2
			while (i<=nrow(mergemat))
	{
		helpvec1<-mergemat[i,][1]
				helpvec2<-mergemat[i,][2]
				if (helpvec1>0)
		{
			helpvec1<-na.exclude(myorder[helpvec1,])
		}
		if (helpvec2>0)
		{
			helpvec2<-na.exclude(myorder[helpvec2,])
		}
		myorder[i,]<-na.add(c(helpvec1,helpvec2))
				i<-i+1
	}
	return(-myorder[i-1,])
}

if (effect.family=="gaussian")
{
if (var.select)
{
if (!(length(transformed.par)==6)) {stop("transformed.par is of a wrong size")}
	{
	btree<-bclustvsG(y,repno,transformed.par)
	}
}
if (!var.select)
{
if (!(length(transformed.par)==5)) {stop("transformed.par is of a wrong size")}
	{
	btree<-bclustG(y,repno,transformed.par)
	}
}
}


if (effect.family=="alaplace")
{
if (var.select)
{
if (!(length(transformed.par)==7)) {stop("transformed.par is of a wrong size")}
	{
	btree<-bclustvsAL(y,repno,transformed.par)
	}
}
if (!var.select)
{
if (!(length(transformed.par)==6)) {stop("transformed.par is of a wrong size")}
	{
	btree<-bclustAL(y,repno,transformed.par)
	}
}
}
if (sum(is.na(btree$height))>0) {stop("dendrogram heights includes NAs, transformed.par values may be inadequately adjusted")}
bclust.tree<-btree
bclust.tree$logposterior<- -btree$height
minheight<-min(btree$height)
minindex<-which(btree$height==minheight)
increment<-diff(btree$height,lag=1)
bclust.tree$height<-diffinv(abs(increment),lag=1)
bclust.tree$clust.number<-((length(repno)-1):1)
bclust.tree$cut<-bclust.tree$height[minindex]
bclust.tree$optim.alloc<-cutree(bclust.tree,h=bclust.tree$cut)
bclust.tree$optim.clustno<-max(bclust.tree$optim.alloc)
bclust.tree$data<-y
bclust.tree$repno<-repno
bclust.tree$transformed.par<-transformed.par
bclust.tree$labels<-labels
bclust.tree$effect.family<-effect.family
bclust.tree$var.select<-var.select

oldClass(bclust.tree)<-c("bclustvs","hclust")
return(bclust.tree)
}



loglikelihood<-function(x.mean,x.css,
repno,transformed.par,effect.family="gaussian",var.select=TRUE)
{
n<-length(as.vector(x.mean))
if (effect.family=="gaussian")
{
	if(var.select)
	{
	if (!(length(transformed.par)==6))  {stop("transformed.par is of a wrong size")}
	return(-.C("RlogmargtmvsG",PACKAGE="bclust",
 	theta=as.double(transformed.par), 
 	mean= as.double(x.mean), 
  	css=as.double(x.css),
  	nrow=as.integer(nrow(x.css)),
   	ncol=as.integer(ncol(x.css)),r=as.double(repno),
    	nr=as.integer(length(repno)), result=as.double(0))$result)
	} else
	{
	if (!(length(transformed.par)==5))  {stop("transformed.par is of a wrong size")}
	return(-.C("RlogmargtmG",PACKAGE="bclust",
 	theta=as.double(transformed.par), 
   	mean= as.double(x.mean),
   	css=as.double(x.css), n=as.integer(n),
   	r=as.double(repno), nr=as.integer(length(repno)),
   	result=as.double(0))$result)
	}
}
if (effect.family=="alaplace")
{
	if(var.select)
	{
	if (!(length(transformed.par)==7))  {stop("transformed.par is of a wrong size")}
	return(-.C("RlogmargtmvsAL",PACKAGE="bclust",
 	theta=as.double(transformed.par), 
  	mean= as.double(as.vector(x.mean)),
  	css=as.double(as.vector(x.css)), 
	nrow=as.integer(nrow(x.css)), ncol=as.integer(ncol(x.css)),
  	r=as.double(repno), nr=as.integer(length(repno)),
  	result=as.double(0))$result)
	} else
	{
	if (!(length(transformed.par)==6))  {stop("transformed.par is of a wrong size")}
	return(-.C("RlogmargtmAL",PACKAGE="bclust",
 	theta=as.double(transformed.par), 
 	mean= as.double(x.mean),
 	css=as.double(x.css), n=as.integer(n), r=as.double(repno),
  	nr=as.integer(length(repno)), result=as.double(0))$result)
	}
 stop("effect.family is  not in the list")	
}
}



meancss<-function(x,rep.id=1:nrow(x))
{
id.order<-order(rep.id)
repno<-as.vector(table(rep.id[id.order]))
y<-x[id.order,]
gendiag<-function(x.vec)
{
y<-rep(sum(x.vec),2*length(x.vec)-1)
y[2*1:length(x.vec)-1]<-x.vec
matrix(rep( rep(1:0,length=2*length(x.vec)-1),y),sum(x.vec),length(x.vec))
}
matsum<-function(x.mat,rowsby) {t(t(x.mat)%*%gendiag(rowsby))}
matrep<-function(x.mat,eachrow)
{
result<-c()
for (i in 1:nrow(x.mat))
result<-c(result,rep(x.mat[i,],eachrow[i]))
return(matrix(result,ncol=ncol(x.mat),byrow=TRUE))
}

R<-matrix(rep(repno,ncol(y)),ncol=ncol(y))
x.mean<-matsum(y,repno)/R
x.css<-matsum((y-matrep(x.mean,repno))^2,repno) 
return(list(mean=x.mean,css=x.css,repno=repno))
}

imp<-function(x)
{
typenomaker<-function(label)
{
helpvec<-c()
for (i in 1:max(label))
      {
      helpvec[i]<-sum(label==i)
      }
return(helpvec)
}

importancevsAL<-function(y,repno,label,hyperparameters)
{
typeno<-typenomaker(label)
lvs1<-matrix(.C("Rlogmarg1datavsAL", PACKAGE="bclust",as.double(y), as.integer(nrow(y)), as.integer(ncol(y)), as.double(repno) ,
as.integer(length(repno)), as.double(typeno), as.integer(length(typeno)), as.double(hyperparameters),result=as.double(rep(0,length(typeno)*ncol(y) )  ) )$result,ncol=ncol(y))
lvs0<-matrix(.C("Rlogmarg0datavsAL", PACKAGE="bclust",as.double(y), as.integer(nrow(y)), as.integer(ncol(y)), as.double(repno) ,
as.integer(length(repno)), as.double(typeno), as.integer(length(typeno)), as.double(hyperparameters),result=as.double(rep(0,length(typeno)*ncol(y) )  ) )$result,ncol=ncol(y))
return(apply(lvs1,2,sum)-apply(lvs0,2,sum) )
}

importanceAL<-function(y,repno,label,hyperparameters)
{
typeno<-typenomaker(label)
lvs1<-matrix(.C("Rlogmarg1dataAL", PACKAGE="bclust",as.double(y), as.integer(nrow(y)), as.integer(ncol(y)), as.double(repno) ,
as.integer(length(repno)), as.double(typeno), as.integer(length(typeno)), as.double(hyperparameters),result=as.double(rep(0,length(typeno)*ncol(y) )  ) )$result,ncol=ncol(y))
lvs0<-matrix(.C("Rlogmarg0dataAL", PACKAGE="bclust",as.double(y), as.integer(nrow(y)), as.integer(ncol(y)), as.double(repno) ,
as.integer(length(repno)), as.double(typeno), as.integer(length(typeno)), as.double(hyperparameters),result=as.double(rep(0,length(typeno)*ncol(y) )  ) )$result,ncol=ncol(y))
return(lvs1-lvs0)
}


importanceG<-function(y,repno,label,hyperparameters)
{
typeno<-typenomaker(label)
lvs1<-matrix(.C("Rlogmarg1dataG", PACKAGE="bclust",as.double(y), as.integer(nrow(y)), as.integer(ncol(y)), as.double(repno) ,
as.integer(length(repno)), as.double(typeno), as.integer(length(typeno)), as.double(hyperparameters),result=as.double(rep(0,length(typeno)*ncol(y) )  ) )$result,ncol=ncol(y))
lvs0<-matrix(.C("Rlogmarg0dataG", PACKAGE="bclust",as.double(y), as.integer(nrow(y)), as.integer(ncol(y)), as.double(repno) ,
as.integer(length(repno)), as.double(typeno), as.integer(length(typeno)), as.double(hyperparameters),result=as.double(rep(0,length(typeno)*ncol(y) )  ) )$result,ncol=ncol(y))
return(lvs1-lvs0)
}

importancevsG<-function(y,repno,label,hyperparameters)
{
typeno<-typenomaker(label)
lvs1<-matrix(.C("Rlogmarg1datavsG", PACKAGE="bclust",as.double(y), as.integer(nrow(y)), as.integer(ncol(y)), as.double(repno) ,
as.integer(length(repno)), as.double(typeno), as.integer(length(typeno)), as.double(hyperparameters),result=as.double(rep(0,length(typeno)*ncol(y) )  ) )$result,ncol=ncol(y))
lvs0<-matrix(.C("Rlogmarg0datavsG", PACKAGE="bclust",as.double(y), as.integer(nrow(y)), as.integer(ncol(y)), as.double(repno) ,
as.integer(length(repno)), as.double(typeno), as.integer(length(typeno)), as.double(hyperparameters),result=as.double(rep(0,length(typeno)*ncol(y) )  ) )$result,ncol=ncol(y))
return(apply(lvs1,2,sum)-apply(lvs0,2,sum) )
}




if (!(class(x)[1]=="bclustvs")) {stop("x should be a bclustvs object")}

optimlab<-cutree(x,h=x$cut)
if (x$optim.clustno<2) {stop("all data in one cluster, importances are useless")}

if (x$var.select) 
{
	if(x$effect.family=="gaussian")
	{
	varimp<-importancevsG(x$data[order(rep(optimlab,x$repno)),],
	x$repno,x$optim.alloc,x$transformed.par)
	} else
	{
	if(x$effect.family=="alaplace")
	{
	varimp<-importancevsAL(x$data[order(rep(x$optim.alloc,x$repno)),],
	x$repno,x$optim.alloc,x$transformed.par)
	}
	}
} else
	{warning("the model is without variable selection, just the variable-cluster importance is
	calculated");varimp<-NULL;varlabel<-NULL;varorder<-NULL}
varorder<-order(varimp,decreasing=TRUE)
#varimp<-sort(varimp,decreasing=TRUE)
if (is.null(colnames(x$data))) {varlabel<-paste(1:ncol(x$data))} else {varlabel<-colnames(x$data)}
#varlabel<-varlabel[varorder]
	if(x$effect.family=="gaussian")
	{
	varclustimp<-importanceG(x$data[order(rep(optimlab,x$repno)),],
	x$repno,x$optim.alloc,x$transformed.par)
	} else
	{
	if(x$effect.family=="alaplace")
	{
	varclustimp<-importanceAL(x$data[order(rep(optimlab,x$repno)),],
	x$repno,x$optim.alloc,x$transformed.par)
	}
	}
#varclustimp<-varclustimp[,varorder]
return(list(var=varimp,varclust=varclustimp,repno=x$repno,labels=varlabel,order=varorder))
}


viplot<-function(varimp,xlab,xlab.mar=5,ylab.mar=4,xlab.srt=90,xlab.cex=1,
sort=FALSE,...)
{
par(mar = c(xlab.mar,ylab.mar, 0.5, 0) + 0.1) # leave some space for labels
if(missing(varimp)) stop ("varimp is missing")
if(missing(xlab)) {xlab<-paste(1:length(varimp))}
if(!(is.vector(xlab)))stop("xlab is not a vector")
if(!(is.vector(varimp)))stop("var is not a vector")
if(!(is.numeric(varimp)))stop("var is not numeric")
if(!((xlab.mar>=0)|(ylab.mar>=0)))stop("margin value is not numeric")
if(!(is.numeric(xlab.srt)))stop("xlab.srt is not numeric")
if(!(xlab.cex>0))stop("xlab.cex is not appropriate")
if(!sort)
{
bp <- barplot(varimp,...) #plot bars
text(bp, par("usr")[3] - 0.5, srt = xlab.srt, adj = 1,
     labels = xlab, xpd = TRUE,cex=xlab.cex,) #plot variable labels
} else {
sortorder<-order(varimp,decreasing=TRUE)
bp <- barplot(varimp[sortorder],...) #plot bars
text(bp, par("usr")[3] - 0.5, srt = xlab.srt, adj = 1,
     labels = xlab[sortorder], xpd = TRUE,cex=xlab.cex) 
	}
}



teethplot<-function(x,teeth.space=0.25,teeth.lwd=1)
{
label.dr<-x$optim.alloc[order.dendrogram(as.dendrogram(x))]
if ((teeth.space>0.25)| (teeth.space<0)) stop("teeth.space is wrongly adjusted")
if (teeth.lwd<0) stop("teeth.lwd is wrongly adjusted")
breakpoints<-function(label)
{
result<-c()
j<-1
for (i in 1:(length(label)-1))
    {
    if (!(label[i]==label[i+1])) {result[j]<-i;j<-j+1}
    }
return(result)
}
typenomaker<-function(label)
{
helpvec<-c()
for (i in 1:max(label))
      {
      helpvec[i]<-sum(label==i)
      }
return(helpvec)
}
relabel<-function(current.label)
{
help.vec<-rep(NA,length(current.label))
j<-1
labelvalue<-current.label[1]
help.vec [which(current.label==labelvalue)]<-j
current.label[which(current.label==labelvalue)]<-0

for (i in 2:length(current.label))
        {
        if ((current.label[i]!=labelvalue) & (current.label[i]!=0))
                {
                labelvalue<-current.label[i]
                j<-j+1
                help.vec [which(current.label==labelvalue)]<-j
                current.label[which(current.label==labelvalue)]<-0
                } 
        }
return(help.vec)
}
label.dr<-relabel(label.dr)
breakpoints.dr<-breakpoints(label.dr)
image(1:3, 1:length(label.dr), matrix(0,3,length(label.dr)), xlim = c(0, 
        1.2 ), ylim = c(0.5,length(label.dr) + 0.5), axes=FALSE,xlab = "", 
        ylab = "", col="white") #this makes plots comparable in axes

points(c(0,1),c(0.5+(0.5-teeth.space),0.5+(0.5-teeth.space)),
type="l",lwd=teeth.lwd)
points(c(0,1),c(length(label.dr)+teeth.space,
length(label.dr)+teeth.space),type="l",lwd=teeth.lwd)

for (i in 1:length(breakpoints.dr))
{
points(c(0,1),c(breakpoints.dr[i]+0.5+teeth.space,breakpoints.dr[i]+0.5+teeth.space),
type="l",lwd=teeth.lwd)
points(c(0,1),c(breakpoints.dr[i]+0.5-teeth.space,breakpoints.dr[i]+0.5-teeth.space),
type="l",lwd=teeth.lwd)
}

typeno.dr<-typenomaker(label.dr)
points(c(1,1),c(sum(typeno.dr)+teeth.space,
sum(typeno.dr)-typeno.dr[length(typeno.dr)]+0.5+teeth.space),type="l",lwd=teeth.lwd)

points(c(1,1),c(0.5+(0.5-teeth.space),
typeno.dr[1]+0.5-teeth.space),type="l",lwd=teeth.lwd)

if (length(typeno.dr)>2)
 {
for (i in 2:(length(typeno.dr)-1))
  {
  points(c(1,1),c(0.5+breakpoints.dr[i-1]+teeth.space,
  breakpoints.dr[i]+0.5-teeth.space),type="l",lwd=teeth.lwd)
  }
 }
}



profileplot<-function(x,rep.id,labels=NULL,scale=1,col.names=colnames(x),plot.order=1:max(rep.id),
xlab.mar=5,ylab.mar=5,xlab.cex=0.8,ylab.cex=1,profile.col=rep(1,max(rep.id)),
blob.matrix=matrix(0,ncol=ncol(x),nrow=max(rep.id)),blob.col=rev(heat.colors(max(blob.matrix))),blob.cex=1.5)
{
if (!(length(plot.order)==max(rep.id))) stop("plot.order does not match with rep.id")
if (missing(x)){stop("data are not specified")}
if (!is.matrix(x)){stop("data must be a matrix")}
if (!(length(rep.id)==nrow(x))){stop("rep.id mismatches")}
if ((sum(is.na(x))>0) | (sum(is.na(rep.id))>0)) stop("NA is not allowed")
if (!(length(blob.col)==max(blob.matrix))) {stop("blob.matrix does not match with blob.col")}
id.order<-order(rep.id)
repno<-as.vector(table(rep.id[id.order]))

if (is.null(labels)){if (is.null(rownames(x))){labels<-paste(1:length(repno))} else  
{if (length(repno)==nrow(x)){labels<-rownames(x)};labels<-labels[id.order]} } else 
{if (!(length(labels)==length(repno))){warning("label length does not match with the data")}}

type.names<-labels

gendiag<-function(x.vec)
{
y<-rep(sum(x.vec),2*length(x.vec)-1)
y[2*1:length(x.vec)-1]<-x.vec
matrix(rep( rep(1:0,length=2*length(x.vec)-1),y),sum(x.vec),length(x.vec))
}
matsum<-function(x.mat,rowsby) {t(t(x.mat)%*%gendiag(rowsby))}


y<-x[order(rep.id),]
ybar<-(matsum(y,repno)/repno)
y<-y/scale
ybar<-ybar/scale
plot.rowfinder<-function(repno,i,j)
{
if (i==1) {return(j)} else{return(sum(repno[1:(i-1)])+j)} 
}
xlab<-col.names
ylab<-type.names
repno.mat<-matrix(rep(repno,ncol(y)),ncol=ncol(y))
ybar<-matsum(y,repno)/repno.mat
#if (is.null(xlab)) {profile.bmargin<-0} else #conflicts inside bclust plot
profile.bmargin<-xlab.mar #bottom margin for xlabels
#if (is.null(ylab)) {profile.rmargin<-0.2} else 
profile.rmargin<-ylab.mar #right margin for ylabels #conflicts inside bclust plot
par(mar=c(profile.bmargin,0,0,profile.rmargin))
image(1:ncol(y), 1:length(repno), t(matrix(0,ncol=ncol(y),nrow=length(repno))), axes = FALSE, xlim = c(0.5, ncol(y) + 0.5), ylim = c(0.5, length(repno) + 0.5), xlab = "", 
       ylab = "", col="white")

for (i in 1:length(plot.order))
{
abline(h=i,col='gray')
}
#write labels
if(!is.null(xlab)) {axis(1, 1:ncol(y), las = 2, line = -0.5, tick = 0, 
            labels = colnames(y), cex.axis = xlab.cex)}
if(!is.null(ylab)){axis(4, 1:length(repno), las = 2, line = -0.5, tick = 0, 
        labels = type.names[plot.order], cex.axis =ylab.cex)} #r.cex)

## Importance shoule add here


for (i in 1:length(plot.order))
{
  for (j in 1:repno[i])
  {
  points(1:ncol(y),i+y[plot.rowfinder(repno,plot.order[i],j),],type='l',col=profile.col[i])
  }
}

for (i in 1:length(plot.order))
{
  for (j in 1:ncol(y))
  {
  if(blob.matrix[plot.order[i],j]>0) points(j,i+ybar[plot.order[i],j],pch=16,cex=blob.cex,col=blob.col[blob.matrix[plot.order[i],j]])
  }
}
}


ditplot<-function(x,xlab=colnames(x$data),ylab=x$labels,
xlab.cex=1,ylab.cex=1,dendrogram.lwd=1,dendrogram.size=2,xlab.mar=3,ylab.mar=3, 
image.col=rainbow(20), horizbar.plot=FALSE,horizbar.col=rev(c(heat.colors(5)[-4],"white")), horizbar.distance=4,varimp=rep(0,ncol(x$data)),horizbar.size=0.5,vertbar=NULL,
vertbar.col=rainbow(max(vertbar)),teeth.size=0.25,plot.width=10) {
relabel<-function(current.label)
{
help.vec<-rep(NA,length(current.label))
j<-1
labelvalue<-current.label[1]
help.vec [which(current.label==labelvalue)]<-j
current.label[which(current.label==labelvalue)]<-0

for (i in 2:length(current.label))
        {
        if ((current.label[i]!=labelvalue) & (current.label[i]!=0))
                {
                labelvalue<-current.label[i]
                j<-j+1
                help.vec [which(current.label==labelvalue)]<-j
                current.label[which(current.label==labelvalue)]<-0
                } 
        }
return(help.vec)
}

gendiag<-function(x.vec)
{
y<-rep(sum(x.vec),2*length(x.vec)-1)
y[2*1:length(x.vec)-1]<-x.vec
matrix(rep( rep(1:0,length=2*length(x.vec)-1),y),sum(x.vec),length(x.vec))
}
matsum<-function(x.mat,rowsby) {t(t(x.mat)%*%gendiag(rowsby))}
logbfcoder<-function(x)
{
result<-c()
for (i in 1:length(x))
   {
           if (x[i]<= 0 ) {result[i]<- 0}
                    else{if (x[i]< 1 ) {result[i]<- 1}
                          else{if (x[i]< 3 ) {result[i]<- 2}
                                else{if (x[i]< 5 ) {result[i]<- 3}
                                      else{result[i]<- 4}
                                    }
                               }     
                        }
                   
     }
return(result)
}

if(is.null(varimp)) {stop("varimp is specified NULL")}

x.data<-(matsum(x$data,x$repno)/x$repno)

cutplot.dendrogramh = function(x, h, cluscol=NULL, leaflab= "none", horiz=TRUE, lwd=1, ...)
{
    if (missing(h)) {
        return(plot(x, leaflab=leaflab, ...))
    }
    
    # Not nice, but necessary
    pn  = plotnode
    
    opar = par()[c("col","lwd")]
    on.exit(par(opar))
    par(lwd=lwd)
   
    x = cut(x, h)
    plot(x[[1]],axes=FALSE, horiz=TRUE,leaflab="none", yaxs="i",xaxs="i",...)
    
    x = x[[2]]
    K = length(x)
    if (is.null(cluscol)) {
       cluscol = rainbow(K)
    }
    x1 = 1
    for (k in 1:K) {
        x2 = x1 + attr(x[[k]],"members")-1
        par(col=cluscol[k])
        pn(x1,x2, x[[k]], type="rectangular", center=FALSE, 
                 leaflab=leaflab, nodePar=NULL, edgePar=list(), horiz=horiz)
        x1 = x2 + 1
   }
abline(v=h,col="gray",lwd=lwd,lty=2)   
}


layout(matrix(c(2,1,3,4,0,5,0,0),2,4,byrow=TRUE), c(dendrogram.size,plot.width,teeth.size,0.12), c(5.3,horizbar.size), respect=TRUE) # defines the space of plot

if (is.null(xlab)) {image.bmargin<-0} else {image.bmargin<-xlab.mar} #bottom margin for xlabels
if (is.null(ylab)) {image.rmargin<-0.2} else {image.rmargin<-ylab.mar} #right margin for ylabels

            bclust.dendro <- as.dendrogram(x)
           rowInd <- order.dendrogram(bclust.dendro)
            varimp.order<-order(varimp,decreasing=TRUE)
            x.data<-x.data[rowInd,varimp.order] 
            ylab<-ylab[rowInd]
            xlab<-xlab[varimp.order]
            varimp<-sort(varimp,decreasing=TRUE)
            varimp<-logbfcoder(varimp)
# plot 1
par(mar=c(horizbar.distance+image.bmargin,0,0,image.rmargin))
image(1:ncol(x.data), 1:nrow(x.data), t(x.data), axes = FALSE, xlim = c(0.5, 
       ncol(x.data) + 0.5), ylim = c(0.5, nrow(x.data) + 0.5), xlab = "", 
       ylab = "", col=image.col)
#writes labels
if(!is.null(xlab)) {axis(1, 1:ncol(x.data), las = 2, line = -0.5, tick = 0, 
            labels = xlab, cex.axis = xlab.cex)}
if(!is.null(ylab)){axis(4, 1:nrow(x.data), las = 2, line = -0.5, tick = 0, 
        labels = ylab, cex.axis =ylab.cex)} #r.cex)

par(mar=c(horizbar.distance+image.bmargin,0.5,0,0))
# plot 2
cutplot.dendrogramh(bclust.dendro, h=x$cut,  horiz=TRUE,leaflab = "none",lwd=dendrogram.lwd)




#plot 3
    par(mar=c(horizbar.distance+image.bmargin,0,0,0))
   label.dr<-relabel(x$optim.alloc[rowInd])
    teethplot(x)

# plot 4
if (!(is.null(vertbar)))
{
par(mar=c(horizbar.distance+image.bmargin,0,0,0))
vertbar<-vertbar[rowInd] 
image(1,1:nrow(x.data),t(matrix(vertbar,nrow(x.data))),axes=FALSE,ylim = c(0.5, nrow(x.data) + 0.5),xlab="",ylab="",col=vertbar.col)
} else {image(1,1:nrow(x.data),t(matrix(0,nrow(x.data))),axes=FALSE,ylim = c(0.5, nrow(x.data) + 0.5),xlab="",ylab="",col="white")}
 #plot 5   
if (horizbar.plot)
{
    par(mar=c(0,0,0,image.rmargin))
    image(1:ncol(x.data),1,matrix(varimp,ncol=1), axes = FALSE, xlab = "",ylab ="",col=horizbar.col)
}
}








dptplot<-function(x,scale=1,xlab=colnames(x$data),ylab=x$labels,
xlab.cex=1,ylab.cex=1,dendrogram.lwd=1,dendrogram.size=2,xlab.mar=3,ylab.mar=3,  horizbar.plot=FALSE,horizbar.col=rev(c(heat.colors(5)[-4],"white")), horizbar.distance=4,varimp=rep(0,ncol(x$data)),horizbar.size=0.5,vertbar=NULL,
vertbar.col=rainbow(max(vertbar)),teeth.size=0.25,plot.width=10) {
relabel<-function(current.label)
{
help.vec<-rep(NA,length(current.label))
j<-1
labelvalue<-current.label[1]
help.vec [which(current.label==labelvalue)]<-j
current.label[which(current.label==labelvalue)]<-0

for (i in 2:length(current.label))
        {
        if ((current.label[i]!=labelvalue) & (current.label[i]!=0))
                {
                labelvalue<-current.label[i]
                j<-j+1
                help.vec [which(current.label==labelvalue)]<-j
                current.label[which(current.label==labelvalue)]<-0
                } 
        }
return(help.vec)
}

gendiag<-function(x.vec)
{
y<-rep(sum(x.vec),2*length(x.vec)-1)
y[2*1:length(x.vec)-1]<-x.vec
matrix(rep( rep(1:0,length=2*length(x.vec)-1),y),sum(x.vec),length(x.vec))
}
matsum<-function(x.mat,rowsby) {t(t(x.mat)%*%gendiag(rowsby))}
logbfcoder<-function(x)
{
result<-c()
for (i in 1:length(x))
   {
           if (x[i]<= 0 ) {result[i]<- 0}
                    else{if (x[i]< 1 ) {result[i]<- 1}
                          else{if (x[i]< 3 ) {result[i]<- 2}
                                else{if (x[i]< 5 ) {result[i]<- 3}
                                      else{result[i]<- 4}
                                    }
                               }     
                        }
                   
     }
return(result)
}


cutplot.dendrogramh = function(x, h, cluscol=NULL, leaflab= "none", horiz=TRUE, lwd=1, ...)
{
    if (missing(h)) {
        return(plot(x, leaflab=leaflab, ...))
    }
    
    # Not nice, but necessary
    pn  = plotnode
    
    opar = par()[c("col","lwd")]
    on.exit(par(opar))
    par(lwd=lwd)
   
    x = cut(x, h)
    plot(x[[1]],axes=FALSE, horiz=TRUE,leaflab="none", yaxs="i",xaxs="i",...)
    
    x = x[[2]]
    K = length(x)
    if (is.null(cluscol)) {
       cluscol = rainbow(K)
    }
    x1 = 1
    for (k in 1:K) {
        x2 = x1 + attr(x[[k]],"members")-1
        par(col=cluscol[k])
        pn(x1,x2, x[[k]], type="rectangular", center=FALSE, 
                 leaflab=leaflab, nodePar=NULL, edgePar=list(), horiz=horiz)
        x1 = x2 + 1
   }
abline(v=h,lwd=lwd,col="gray",lty=2)   
}


layout(matrix(c(2,1,3,4,0,5,0,0),2,4,byrow=TRUE), c(dendrogram.size,plot.width,teeth.size,0.12), c(5.3,horizbar.size), respect=TRUE) # defines the space of plot

if (is.null(xlab)) {image.bmargin<-0} else {image.bmargin<-xlab.mar} #bottom margin for xlabels
if (is.null(ylab)) {image.rmargin<-0.2} else {image.rmargin<-ylab.mar} #right margin for ylabels

            bclust.dendro <- as.dendrogram(x)
           rowInd <- order.dendrogram(bclust.dendro)
           varimp.order<-order(varimp,decreasing=TRUE)

            x.data<-x$data[,varimp.order]
            xlab<-xlab[varimp.order]
            varimp<-sort(varimp,decreasing=TRUE)
            varimp<-logbfcoder(varimp)
# plot 1

profileplot(x=x.data, rep.id=rep(1:length(x$repno),x$repno), labels = ylab, scale = scale, col.names = xlab, plot.order = rowInd, xlab.mar =(horizbar.distance+image.bmargin), ylab.mar = image.rmargin, xlab.cex = xlab.cex, ylab.cex = ylab.cex)


par(mar=c(horizbar.distance+image.bmargin,0.5,0,0))
# plot 2
cutplot.dendrogramh(bclust.dendro, h=x$cut,  horiz=TRUE,leaflab = "none",lwd=dendrogram.lwd)




#plot 3
    par(mar=c(horizbar.distance+image.bmargin,0,0,0))
   label.dr<-relabel(x$optim.alloc[rowInd])
    teethplot(x)

# plot 4
if (!(is.null(vertbar)))
{
par(mar=c(horizbar.distance+image.bmargin,0,0,0))
vertbar<-vertbar[rowInd] 
image(1,1:nrow(x.data),t(matrix(vertbar,nrow(x.data))),axes=FALSE,ylim = c(0.5, nrow(x.data) + 0.5),xlab="",ylab="",col=vertbar.col)
} else {image(1,1:nrow(x.data),t(matrix(0,nrow(x.data))),axes=FALSE,ylim = c(0.5, nrow(x.data) + 0.5),xlab="",ylab="",col="white")}
 #plot 5   
if (horizbar.plot)
{
    par(mar=c(0,0,0,image.rmargin))
    image(1:ncol(x.data),1,matrix(varimp,ncol=1), axes = FALSE, xlab = "",ylab ="",col=horizbar.col)
}
}




















bdiscrim<-function(training, training.id=NULL, training.labels=NULL, predict=NULL, predict.label=rownames(predict)[1],  effect.family="gaussian",var.select=TRUE,
transformed.par,priorprob=rep(1,max(training.id)+1))
{

typenomaker<-function(label)
{
helpvec<-c()
for (i in 1:max(label))
      {
      helpvec[i]<-sum(label==i)
      }
return(helpvec)
}

importancevsAL<-function(y,repno,label,hyperparameters)
{
typeno<-typenomaker(label)
lvs1<-matrix(.C("Rlogmarg1datavsAL", PACKAGE="bclust",as.double(y), as.integer(nrow(y)), as.integer(ncol(y)), as.double(repno) ,
as.integer(length(repno)), as.double(typeno), as.integer(length(typeno)), as.double(hyperparameters),result=as.double(rep(0,length(typeno)*ncol(y) )  ) )$result,ncol=ncol(y))
lvs0<-matrix(.C("Rlogmarg0datavsAL", PACKAGE="bclust",as.double(y), as.integer(nrow(y)), as.integer(ncol(y)), as.double(repno) ,
as.integer(length(repno)), as.double(typeno), as.integer(length(typeno)), as.double(hyperparameters),result=as.double(rep(0,length(typeno)*ncol(y) )  ) )$result,ncol=ncol(y))
return(apply(lvs1,2,sum)-apply(lvs0,2,sum) )
}

importanceAL<-function(y,repno,label,hyperparameters)
{
typeno<-typenomaker(label)
lvs1<-matrix(.C("Rlogmarg1dataAL", PACKAGE="bclust",as.double(y), as.integer(nrow(y)), as.integer(ncol(y)), as.double(repno) ,
as.integer(length(repno)), as.double(typeno), as.integer(length(typeno)), as.double(hyperparameters),result=as.double(rep(0,length(typeno)*ncol(y) )  ) )$result,ncol=ncol(y))
lvs0<-matrix(.C("Rlogmarg0dataAL", PACKAGE="bclust",as.double(y), as.integer(nrow(y)), as.integer(ncol(y)), as.double(repno) ,
as.integer(length(repno)), as.double(typeno), as.integer(length(typeno)), as.double(hyperparameters),result=as.double(rep(0,length(typeno)*ncol(y) )  ) )$result,ncol=ncol(y))
return(lvs1-lvs0)
}


importanceG<-function(y,repno,label,hyperparameters)
{
typeno<-typenomaker(label)
lvs1<-matrix(.C("Rlogmarg1dataG", PACKAGE="bclust",as.double(y), as.integer(nrow(y)), as.integer(ncol(y)), as.double(repno) ,
as.integer(length(repno)), as.double(typeno), as.integer(length(typeno)), as.double(hyperparameters),result=as.double(rep(0,length(typeno)*ncol(y) )  ) )$result,ncol=ncol(y))
lvs0<-matrix(.C("Rlogmarg0dataG", PACKAGE="bclust",as.double(y), as.integer(nrow(y)), as.integer(ncol(y)), as.double(repno) ,
as.integer(length(repno)), as.double(typeno), as.integer(length(typeno)), as.double(hyperparameters),result=as.double(rep(0,length(typeno)*ncol(y) )  ) )$result,ncol=ncol(y))
return(lvs1-lvs0)
}

importancevsG<-function(y,repno,label,hyperparameters)
{
typeno<-typenomaker(label)
lvs1<-matrix(.C("Rlogmarg1datavsG", PACKAGE="bclust",as.double(y), as.integer(nrow(y)), as.integer(ncol(y)), as.double(repno) ,
as.integer(length(repno)), as.double(typeno), as.integer(length(typeno)), as.double(hyperparameters),result=as.double(rep(0,length(typeno)*ncol(y) )  ) )$result,ncol=ncol(y))
lvs0<-matrix(.C("Rlogmarg0datavsG", PACKAGE="bclust",as.double(y), as.integer(nrow(y)), as.integer(ncol(y)), as.double(repno) ,
as.integer(length(repno)), as.double(typeno), as.integer(length(typeno)), as.double(hyperparameters),result=as.double(rep(0,length(typeno)*ncol(y) )  ) )$result,ncol=ncol(y))
return(apply(lvs1,2,sum)-apply(lvs0,2,sum) )
}

logclassprobG<-function(ynew,yknown,repnoknown,hyperparameters,logpriorprob)
{
repno<-c(repnoknown,nrow(ynew))
y<-rbind(yknown,ynew)
logclassprob<-c()
for (i in 1:length(repno))
     {
     label<-c(1:(length(repno)-1),i)
     logclassprob[i]<-.C("RloglikbylabelGunif", PACKAGE="bclust",as.double(as.matrix(y)), nrowy=as.integer(nrow(y)), ncoly=as.integer(ncol(y)), 
     repno=as.double(repno), nrepno=as.integer(length(repno)),label=as.double(label),nlabel=as.integer(length(label)),
     theta=as.double(hyperparameters), result=as.double(0))$result+logpriorprob[i]
     }
return(logclassprob)
}


logclassprobvsG<-function(ynew,yknown,repnoknown,hyperparameters,logpriorprob)
{
repno<-c(repnoknown,nrow(ynew))
y<-rbind(yknown,ynew)
logclassprob<-c()
for (i in 1:length(repno))
     {
     label<-c(1:(length(repno)-1),i)
     logclassprob[i]<-.C("RloglikbylabelvsGunif", PACKAGE="bclust",as.double(as.matrix(y)), nrowy=as.integer(nrow(y)), ncoly=as.integer(ncol(y)), 
     repno=as.double(repno), nrepno=as.integer(length(repno)),label=as.double(label),nlabel=as.integer(length(label)),
     theta=as.double(hyperparameters), result=as.double(0))$result+logpriorprob[i]
     }
return(logclassprob)
}


logclassprobAL<-function(ynew,yknown,repnoknown,hyperparameters,logpriorprob)
{
repno<-c(repnoknown,nrow(ynew))
y<-rbind(yknown,ynew)
logclassprob<-c()
for (i in 1:length(repno))
     {
     label<-c(1:(length(repno)-1),i)
     logclassprob[i]<-.C("RloglikbylabelALunif", PACKAGE="bclust",as.double(as.matrix(y)), nrowy=as.integer(nrow(y)), ncoly=as.integer(ncol(y)), 
     repno=as.double(repno), nrepno=as.integer(length(repno)),label=as.double(label),nlabel=as.integer(length(label)),
     theta=as.double(hyperparameters), result=as.double(0))$result+logpriorprob[i]
     }
return(logclassprob)
}

logclassprobvsAL<-function(ynew,yknown,repnoknown,hyperparameters,logpriorprob)
{
repno<-c(repnoknown,nrow(ynew))
y<-rbind(yknown,ynew)
logclassprob<-c()
for (i in 1:length(repno))
     {
     label<-c(1:(length(repno)-1),i)
     logclassprob[i]<-.C("RloglikbylabelvsALunif", PACKAGE="bclust",as.double(as.matrix(y)), nrowy=as.integer(nrow(y)), ncoly=as.integer(ncol(y)), 
     repno=as.double(repno), nrepno=as.integer(length(repno)),label=as.double(label),nlabel=as.integer(length(label)),
     theta=as.double(hyperparameters), result=as.double(0))$result+logpriorprob[i]
     }
return(logclassprob)
}


if (is.null(training.id)){stop("training.id is missing")}
if (missing(training)){stop("training data are not specified")}
if (!is.matrix(training)){stop("data must be a matrix")}
if (missing(transformed.par)){stop("transformed.par is missing")}
if (!is.vector(transformed.par)){stop("transformed.par must be a vector")}
if (!(length(training.id)==nrow(training))){stop("training.id mismatches")}
if ((sum(is.na(training))>0) | (sum(is.na(training.id))>0) | 
(sum(is.na(transformed.par))>0)) stop("NA is not allowed")
if(!(effect.family%in%c("gaussian","alaplace"))) {stop(
"in effect.family just the gaussian and the alaplace families are allowed")}
id.order<-order(training.id)
repno<-as.vector(table(training.id[id.order]))
y<-training[id.order,]
labels<-training.labels
if (is.null(labels)){if (is.null(rownames(training))) {labels<-paste(1:length(repno))} else  {if (length(repno)==nrow(training)){labels<-rownames(training)};labels<-labels[id.order]} } else {if (!(length(labels)==length(repno))){warning("label length does not match with the data")}}

if(!is.null(predict)){
	if (!is.matrix(predict)) {stop("predict should be a matrix")}
	if (!(ncol(predict)==ncol(training))) stop("dimension of predict unmatches with training")
	if(sum(is.na(training))>0) {stop("NA is not allowed in predict")}
	}

bprob<-NULL

if (effect.family=="gaussian")
{
if (var.select)
{
if (!(length(transformed.par)==6)) {stop("transformed.par is of a wrong size")}
	{
	if(!is.null(predict))
		{
		blogprob<-logclassprobvsG(predict,y, repno,transformed.par,log(priorprob))
		blogprob<-matrix(blogprob,nrow=1)
		colnames(blogprob)<-c(labels,"New")
		rownames(blogprob)<-predict.label
		blogprob<-blogprob-max(blogprob,na.rm=TRUE)
		bprob<-exp(blogprob)/(sum(exp(blogprob),na.rm=TRUE))
		} else
		{warning("predict is not specified, just importances are calculated")}
	varimp<-importancevsG(training,repno,1:length(repno),transformed.par)
	varclassimp<-importanceG(training,repno,1:length(repno),transformed.par)
	}
}
if (!var.select)
{
if (!(length(transformed.par)==5)) {stop("transformed.par is of a wrong size")}
	{
	if(!is.null(predict))
		{
		blogprob<-logclassprobG(predict,y, repno,transformed.par,log(priorprob))
		blogprob<-matrix(blogprob,nrow=1)
		colnames(blogprob)<-c(labels,"New")
		rownames(blogprob)<-predict.label
		blogprob<-blogprob-max(blogprob,na.rm=TRUE)
		bprob<-exp(blogprob)/(sum(exp(blogprob),na.rm=TRUE))

		} else
		{warning("predict is not specified, just importances are calculated")}
	warning("no variable selection, just the variable-class importance is calculated")
	varimp<-NULL
	varclassimp<-importanceG(training,repno,1:length(repno),transformed.par)
	}
}
}


if (effect.family=="alaplace")
{
if (var.select)
{
if (!(length(transformed.par)==7)) {stop("transformed.par is of a wrong size")}
	{
	if(!is.null(predict))
		{
		blogprob<-logclassprobvsAL(predict,y, repno,transformed.par,log(priorprob))
		blogprob<-matrix(blogprob,nrow=1)
		colnames(blogprob)<-c(labels,"New")
		rownames(blogprob)<-predict.label
		blogprob<-blogprob-max(blogprob,na.rm=TRUE)
		bprob<-exp(blogprob)/(sum(exp(blogprob),na.rm=TRUE))
		} else
		{warning("predict is not specified, just importances are calculated")}
	varimp<-importancevsAL(training,repno,1:length(repno),transformed.par)
	varclassimp<-importanceAL(training,repno,1:length(repno),transformed.par)
	}
}
if (!var.select)
{
if (!(length(transformed.par)==6)) {stop("transformed.par is of a wrong size")}
	{
	if(!is.null(predict))
		{
		blogprob<-logclassprobAL(predict,y, repno,transformed.par,log(priorprob))
		blogprob<-matrix(blogprob,nrow=1)
		colnames(blogprob)<-c(labels,"New")
		rownames(blogprob)<-predict.label
		blogprob<-blogprob-max(blogprob,na.rm=TRUE)
		bprob<-exp(blogprob)/(sum(exp(blogprob),na.rm=TRUE))
		} else
		{warning("predict is not specified, just importances are calculated")}
	warning("no variable selection, just the variable-class importance is calculated")
	varimp<-NULL
	varclassimp<-importanceAL(training,repno,1:length(repno),transformed.par)
	}
}
}


return(list(probs=bprob,var=varimp,varclass=varclassimp))
}



