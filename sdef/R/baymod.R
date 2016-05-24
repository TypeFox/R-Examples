baymod <-
function(output.ratio,iter=1000,dir=getwd(),conf=95){
load(paste(output.ratio$dataname,".Rdata"))
if(output.ratio$pvalue==FALSE){
data=1-data;
}

#Number of lists
lists = dim(data)[2]

#Calculate all the needed quantities
dim1=dim(data)[1]
ID=seq(1,dim1)
l=length(output.ratio$Common)

#Design
design = designMatrix(lists)
rowdes=nrow(design)-1
row.names(design)<-seq(1,dim(design)[1])
threshold = output.ratio$h
O = matrix(0,l,(rowdes))

for(i in 1:l){
temp = matrix(0,dim1,lists)
for(j in 1:lists){
for(r in 1:dim1){
if(data[r,j]<=threshold[i]) {temp[r,j] = 1}
}
}
xx=designCount(temp,design)
O[i,] <- xx
}

#Dirichlet prior
prior.p=0.001
post.p=matrix(NA,l,(rowdes))

p=array(NA,dim=c(l,(rowdes),iter))
p.s = p
marginal.p=array(NA,dim=c(l,lists,iter))
ratio = matrix(NA,l,iter)

for(i in 1:l){
for(j in 1:(rowdes)){
post.p[i,j] <- O[i,j] + prior.p

for(k in 1:iter){   
p.s[i,j,k] <- rgamma(1,post.p[i,j],1)
p[i,j,k] <- p.s[i,j,k]/dim1
}
}
for(j in 1:lists){
for(k in 1:iter){
marginal.p[i,j,k]<-sum(p[i,(as.numeric(row.names(design[design[,j]==1,]))-1),k])
}
for(k in 1:iter){
ratio[i,k]<-p[i,rowdes,k]/prod(marginal.p[i,,k])
}
}
}
#########################################
#CI for ratio
quantile = matrix(NA,l,3)
lower=(100-conf)/200
upper=1-lower
for(i in 1:l){
    quantile[i,1] <- quantile(ratio[i,],lower,na.rm=TRUE)
    quantile[i,2] <- quantile(ratio[i,],0.5,na.rm=TRUE)
    quantile[i,3] <- quantile(ratio[i,],upper,na.rm=TRUE)
}

dimnames(quantile)[[2]]<-c(paste(as.character(lower*100),"%"),"Median",paste(as.character(upper*100),"%"))

lim1<-matrix(0,l,2)
for (i in 1:l) {
lim1[i,1]<-quantile[i,1]
lim1[i,2]<-quantile[i,3]}
y1<-seq(1:l)
y1<-matrix(y1,l,2)
if(length(quantile[round(quantile[,1],3)>1,2])>0){
Rmax = max(quantile[round(quantile[,1],3)>1,2])
hmax = output.ratio$h[quantile[,2]==Rmax]
}
if(length(quantile[round(quantile[,1],3)>1,2])==0){
Rmax=1
}

main=paste("Distribution of R(h) with ",as.character(conf),"% credibility interval")

########
par(omd=c(0.1,0.9,0,1))
plot(y1,lim1,xlab="P value",ylab="R",main=main,cex.main=0.7,pch="_",axes=TRUE,yaxt="n",xaxt="n",
ylim=c(0,(max(quantile[1:l,3],na.rm=TRUE)+1*sd(quantile[1:l,3],na.rm=TRUE))),lwd=0.2)
mtext("number of common genes",side=4,line=2,adj=0.5)

if(Rmax==1){
for (i in 1:l) lines(y1[i,],lim1[i,], lty=3,lwd=1.7)
        axis(2, at = c(0,0.5,1,1.5), labels = c(0,0.5,1,1.5), tick = TRUE,cex.axis=0.9)
axis(1, at = seq(1:l), labels = seq(1:l),tick=TRUE,cex=0.9,las=2)
cat("WARNING: the requested contrast is under-represented in the data (Rmax<1)\n")
}

if(Rmax>1){
    if(length(output.ratio$h[quantile[round(quantile[,1],3)>1,2]>=2])>0){
        h2 = max(output.ratio$h[quantile[,2]>=2])
        R2 = quantile[output.ratio$h==h2,2]
if(h2==hmax){
        for (i in 1:l) lines(y1[i,],lim1[i,], lty=3,lwd=1.7)
        axis(2, at = c(0,0.5,1,1.5,Rmax), labels = c(0,0.5,1,1.5,expression(R[max])), tick = TRUE,cex.axis=0.9)
        if(output.ratio$pvalue==TRUE){   
        axis(1, at = c((hmax*100),seq(((hmax*100)+10),100,20)), labels = c(expression(h[max]),seq((hmax+0.1),1,0.2)),tick=TRUE,cex=0.9,las=2)
        }
                if(output.ratio$pvalue==FALSE){   
                axis(1, at = c((hmax*100),seq(((hmax*100)+10),100,20)), labels = c(expression(h[max]),1-seq((h2+0.1),1,0.2)),tick=TRUE,cex=0.9,las=2)
                }
                axis(4, at = c(1,Rmax),labels = c(dim1,output.ratio$Common[output.ratio$h==hmax]),tick=TRUE,cex=0.9)
}
if(h2!=hmax){
        for (i in 1:l) lines(y1[i,],lim1[i,], lty=3,lwd=1.7)
        axis(2, at = c(0,0.5,1,1.5,R2,Rmax), labels = c(0,0.5,1,1.5,expression(R[2]),expression(R[max])), tick = TRUE,cex.axis=0.9)
        if(output.ratio$pvalue==TRUE){   
        axis(1, at = c((hmax*100),(h2*100),seq(((h2*100)+10),100,20)), labels = c(expression(h[max]),expression(h[2]),seq((h2+0.1),1,0.2)),tick=TRUE,cex=0.9,las=2)
        }
                if(output.ratio$pvalue==FALSE){   
                axis(1, at = c((hmax*100),(h2*100),seq(((h2*100)+10),100,20)), labels = c(expression(h[max]),expression(h[2]),1-seq((h2+0.1),1,0.2)),tick=TRUE,cex=0.9,las=2)
                }
                axis(4, at = c(1,R2,Rmax),labels = c(dim1,output.ratio$Common[output.ratio$h==h2],output.ratio$Common[output.ratio$h==hmax]),tick=TRUE,cex=0.9)
        }
abline(h=R2,lty=3,cex=0.7)
}
    if(length(output.ratio$h[quantile[round(quantile[,1],3)>1,2]>=2])==0){

        for (i in 1:l) lines(y1[i,],lim1[i,], lty=3,lwd=1.7)
        axis(2, at = c(seq(0,(Rmax-0.5),0.5),Rmax), labels = c(seq(0,(Rmax-0.5),0.5),expression(R[max])), tick = TRUE,cex.axis=0.9)
        if(output.ratio$pvalue==TRUE){   
        axis(1, at = c((hmax*100),seq(((hmax*100)+10),100,20)), labels = c(expression(h[max]),seq((hmax+0.1),1,0.2)), tick=TRUE,cex=0.9,las=2)
        }
                if(output.ratio$pvalue==FALSE){   
                axis(1, at = c((hmax*100),seq(((hmax*100)+10),100,20)), labels = c(expression(h[max]),1-seq((hmax+0.1),1,0.2)), tick=TRUE,cex=0.9,las=2)
                }
        }

axis(4, at = c(1,Rmax),labels = c(dim1,output.ratio$Common[output.ratio$h==hmax]),tick=TRUE,cex=0.9)
}
abline(h=1,col="black", lwd=1.5)
points(y1[,1],quantile[1:l,2],col="red",cex=0.5)
abline(h=Rmax,lty=3,cex=0.7)
########

ps.options(horizontal=FALSE)
setwd(dir)
postscript("Rh.ps")
par(omd=c(0.1,0.9,0,1))
plot(y1,lim1,xlab="P value",ylab="R",main=main,cex.main=0.7,pch="_",axes=TRUE,yaxt="n",xaxt="n",
ylim=c(0,(max(quantile[1:l,3],na.rm=TRUE)+1*sd(quantile[1:l,3],na.rm=TRUE))),lwd=0.2)
mtext("number of common genes",side=4,line=2,adj=0.5)

if(Rmax==1){
for (i in 1:l) lines(y1[i,],lim1[i,], lty=3,lwd=1.7)
        axis(2, at = c(0,0.5,1,1.5), labels = c(0,0.5,1,1.5), tick = TRUE,cex.axis=0.9)
axis(1, at = seq(1:l), labels = seq(1:l),tick=TRUE,cex=0.9,las=2)
cat("WARNING: the requested contrast is under-represented in the data (Rmax<1)\n")
}

if(Rmax>1){
    if(length(output.ratio$h[quantile[round(quantile[,1],3)>1,2]>=2])>0){
        h2 = max(output.ratio$h[quantile[,2]>=2])
        R2 = quantile[output.ratio$h==h2,2]
if(h2==hmax){
        for (i in 1:l) lines(y1[i,],lim1[i,], lty=3,lwd=1.7)
        axis(2, at = c(0,0.5,1,1.5,Rmax), labels = c(0,0.5,1,1.5,expression(R[max])), tick = TRUE,cex.axis=0.9)
        if(output.ratio$pvalue==TRUE){   
        axis(1, at = c((hmax*100),seq(((hmax*100)+10),100,20)), labels = c(expression(h[max]),seq((hmax+0.1),1,0.2)),tick=TRUE,cex=0.9,las=2)
        }
                if(output.ratio$pvalue==FALSE){   
                axis(1, at = c((hmax*100),seq(((hmax*100)+10),100,20)), labels = c(expression(h[max]),1-seq((h2+0.1),1,0.2)),tick=TRUE,cex=0.9,las=2)
                }
                axis(4, at = c(1,Rmax),labels = c(dim1,output.ratio$Common[output.ratio$h==hmax]),tick=TRUE,cex=0.9)
}
if(h2!=hmax){
        for (i in 1:l) lines(y1[i,],lim1[i,], lty=3,lwd=1.7)
        axis(2, at = c(0,0.5,1,1.5,R2,Rmax), labels = c(0,0.5,1,1.5,expression(R[2]),expression(R[max])), tick = TRUE,cex.axis=0.9)
        if(output.ratio$pvalue==TRUE){   
        axis(1, at = c((hmax*100),(h2*100),seq(((h2*100)+10),100,20)), labels = c(expression(h[max]),expression(h[2]),seq((h2+0.1),1,0.2)),tick=TRUE,cex=0.9,las=2)
        }
                if(output.ratio$pvalue==FALSE){   
                axis(1, at = c((hmax*100),(h2*100),seq(((h2*100)+10),100,20)), labels = c(expression(h[max]),expression(h[2]),1-seq((h2+0.1),1,0.2)),tick=TRUE,cex=0.9,las=2)
                }
                axis(4, at = c(1,R2,Rmax),labels = c(dim1,output.ratio$Common[output.ratio$h==h2],output.ratio$Common[output.ratio$h==hmax]),tick=TRUE,cex=0.9)
        }
abline(h=R2,lty=3,cex=0.7)
}
    if(length(output.ratio$h[quantile[round(quantile[,1],3)>1,2]>=2])==0){

        for (i in 1:l) lines(y1[i,],lim1[i,], lty=3,lwd=1.7)
        axis(2, at = c(seq(0,(Rmax-0.5),0.5),Rmax), labels = c(seq(0,(Rmax-0.5),0.5),expression(R[max])), tick = TRUE,cex.axis=0.9)
        if(output.ratio$pvalue==TRUE){   
        axis(1, at = c((hmax*100),seq(((hmax*100)+10),100,20)), labels = c(expression(h[max]),seq((hmax+0.1),1,0.2)), tick=TRUE,cex=0.9,las=2)
        }
                if(output.ratio$pvalue==FALSE){   
                axis(1, at = c((hmax*100),seq(((hmax*100)+10),100,20)), labels = c(expression(h[max]),1-seq((hmax+0.1),1,0.2)), tick=TRUE,cex=0.9,las=2)
                }
        }

axis(4, at = c(1,Rmax),labels = c(dim1,output.ratio$Common[output.ratio$h==hmax]),tick=TRUE,cex=0.9)
}
abline(h=1,col="black", lwd=1.5)
points(y1[,1],quantile[1:l,2],col="red",cex=0.5)
abline(h=Rmax,lty=3,cex=0.7)
dev.off()

rownames(quantile)<-output.ratio$h

return(quantile=quantile)
}

