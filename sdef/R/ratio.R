ratio <-
function(data,pvalue=TRUE,interval=0.01,name="Distribution of T(h)",dir=getwd(),dataname="dataratio"){

#Define how many lists for the comparison
lists = ncol(data)
dim1=dim(data)[1]

if(pvalue==FALSE){
data=1-data;
}

ID=seq(1,dim1)
threshold = seq(0,1,interval)
l=length(threshold)
int=c()
L=matrix(0,l,lists)

for(i in 1:l){
temp = data<=threshold[i]
for(j in 1:lists){
L[i,j] <- sum(temp[,j])
temp[temp[,j]==FALSE,j]<-0
temp[temp[,j]==TRUE,j]<-1
        }
        int[i] <- sum(apply(temp,1,sum)==lists)
}

#Calculate the ratio for the number of observed features/number of expected ones
expected = apply(L,1,prod)/(dim1)^(lists-1)
ratios = matrix(0,l,1)

for(i in 1:l){
ratios[i,1] <- int[i]/expected[i]
}

#Plot

ratios=ratios[int>0]
thresh.ratios=threshold[int>0]
L = L[int>0,] 
int=int[int>0]

Tmax = max(ratios)
hmax = thresh.ratios[ratios==Tmax]
if(length(hmax)>1){hmax=hmax[1]}
if(Tmax<1){cat("WARNING: the requested contrast is under-represented in the data (Tmax<1)\n")}
ps.options(paper="a4",horizontal=TRUE)
setwd(dir)
ps.options(horizontal=FALSE)
postscript("Th.ps")

main="Distribution of T(h)"
   
par(omd=c(0.1,0.9,0,1))
plot(thresh.ratios,ratios,type="l",
ylab= "T",xlab="P-value",main=main,yaxt="n",xaxt="n",cex.main=0.7,cex.axis=1.2,ylim=c(0,(max(ratios,na.rm=TRUE)+sd(ratios,na.rm=TRUE))))
mtext("number of common genes",side=4,line=2,adj=0.5)

if(Tmax<1.1){
axis(2, at = c(0,0.5,Tmax), labels = c(0,0.5,expression(T[max])), tick = TRUE,cex.axis=0.9)
}
if(Tmax>1.1){
axis(2, at = c(0,0.5,1,Tmax), labels = c(0,0.5,1,expression(T[max])), tick = TRUE,cex.axis=0.9)
}
if(hmax>0.1 & hmax<0.9){
axis(1, at = c(0,hmax,seq((hmax+0.1),1,0.2),1), labels = c(0,expression(h[max]),seq((hmax+0.1),1,0.2),1), tick=TRUE,cex=0.9,las=2)
}
if(hmax<0.1){
axis(1, at = c(hmax,seq((hmax+0.1),1,0.2),1), labels = c(expression(h[max]),seq((hmax+0.1),1,0.2),1), tick=TRUE,cex=0.9,las=2)
}
if(hmax>0.9 & hmax != 1){
axis(1, at = c(0,seq(0.1,(hmax-0.1),0.2),hmax), labels = c(0,(hmax-0.1),expression(h[max])), tick=TRUE,cex=0.9,las=2)
}
if(hmax==1){
axis(1, at = c(seq(0,1,0.2)), labels = c(seq(0,(hmax-0.1),0.2),paste(expression(h[max]),"=1")), tick=TRUE,cex=0.9,las=2)
}

axis(4, at = c(1,Tmax),labels = c(dim1,int[thresh.ratios==hmax]),tick=TRUE,cex=0.9)
abline(h=Tmax,lty=3,cex=0.7)
dev.off() 

save(data,file = paste(dataname,".Rdata"))

########
    
par(omd=c(0.1,0.9,0,1))
plot(thresh.ratios,ratios,type="l",
ylab= "T",xlab="P-value",main=main,yaxt="n",xaxt="n",cex.main=0.7,cex.axis=1.2,ylim=c(0,(max(ratios,na.rm=TRUE)+sd(ratios,na.rm=TRUE))))
mtext("number of common genes",side=4,line=2,adj=0.5)

if(Tmax<1.1){
axis(2, at = c(0,0.5,Tmax), labels = c(0,0.5,expression(T[max])), tick = TRUE,cex.axis=0.9)
}
if(Tmax>1.1){
axis(2, at = c(0,0.5,1,Tmax), labels = c(0,0.5,1,expression(T[max])), tick = TRUE,cex.axis=0.9)
}
if(hmax>0.1 & hmax<0.9){
axis(1, at = c(0,hmax,seq((hmax+0.1),1,0.2),1), labels = c(0,expression(h[max]),seq((hmax+0.1),1,0.2),1), tick=TRUE,cex=0.9,las=2)
}
if(hmax<0.1){
axis(1, at = c(hmax,seq((hmax+0.1),1,0.2),1), labels = c(expression(h[max]),seq((hmax+0.1),1,0.2),1), tick=TRUE,cex=0.9,las=2)
}
if(hmax>0.9 & hmax != 1){
axis(1, at = c(0,seq(0.1,(hmax-0.1),0.2),hmax), labels = c(0,(hmax-0.1),expression(h[max])), tick=TRUE,cex=0.9,las=2)
}
if(hmax==1){
axis(1, at = c(seq(0,1,0.2)), labels = c(seq(0,(hmax-0.1),0.2),paste(expression(h[max]),"=1")), tick=TRUE,cex=0.9,las=2)
}

axis(4, at = c(1,Tmax),labels = c(dim1,int[thresh.ratios==hmax]),tick=TRUE,cex=0.9)
abline(h=Tmax,lty=3,cex=0.7) 

########
col.name=NULL
for(i in 1:lists){
temp=paste("list",i)
col.name=c(col.name,temp)
}

colnames(L)<-(col.name)
rownames(L)<-(thresh.ratios)
ratios=as.matrix(ratios)
rownames(ratios)<-(thresh.ratios)
colnames(ratios)<-"ratio"
int=as.matrix(int)
rownames(int)<-(thresh.ratios)
colnames(int)<-"Genes in common"

return(list(h=thresh.ratios,DE = L, ratios=ratios,Common=int,interval=interval,name=name,pvalue=pvalue,dataname=dataname))

}

