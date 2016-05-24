plotTrait <-
function(dat, traits, focal=dat[,2]){
if (is.vector(dat[,traits]) ==  TRUE) if(is.numeric(dat[,traits]) != TRUE) stop("Trait data must be numeric")
if (is.vector(dat[,traits]) !=  TRUE) for (i in 1:length(dat[,traits])){
if(is.numeric(dat[,traits[i]]) != TRUE) stop("Trait data must be numeric")
}
 
if(sort(unique(focal))[1]!=0) stop('Focals must be a binary string denoting non-focals as 0 and focals as 1')
if(sort(unique(focal))[2]!=1) stop('Focals must be a binary string denoting non-focals as 0 and focals as 1')
if(length(unique(focal))!=2) stop('Focals must be a binary string denoting non-focals as 0 and focals as 1')
 

if(length(traits)>3) stop('Graphical output only available for up to three traits')

if(length(traits)==1){
y<-rep(0,dim(dat)[1])
plot(dat[,traits],y,yaxt='n',xlab=paste(traits),ylab='')
points(dat[which(focal==1),traits],y[which(focal==1)],col=2)
legend('top',horiz=TRUE,inset=c(0,-0.1),xpd=TRUE,bty = "n",c('Non-focal','Focal'),pch=1,col=c(1,2))
}

if(length(traits)==2){
tr1<-traits[1]
tr2<-traits[2]
plot(dat[,tr2],dat[,tr1],ylab=paste(tr1),xlab=paste(tr2)) 
points(dat[which(focal==1),tr1]~dat[which(focal==1),tr2],col=2)
legend('top',horiz=TRUE,inset=c(0,-0.1),xpd=TRUE,bty = "n",c('Non-focal','Focal'),pch=1,col=c(1,2))
}

if(length(traits)==3){
tr1<-traits[1]
tr2<-traits[2]
tr3<-traits[3]
myscat<-scatterplot3d(dat[,tr1],dat[,tr2],dat[,tr3],xlab=paste(tr1),ylab=paste(tr2),zlab=paste(tr3))
myscat$points3d(dat[which(focal==1),tr1],dat[which(focal==1),tr2],dat[which(focal==1),tr3],col=2)
legend('top',horiz=TRUE,inset=c(0,-0.1),xpd=TRUE,bty = "n",c('Non-focal','Focal'),pch=1,col=c(1,2))
}
}
