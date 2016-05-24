#Graphics output of Cliffs delta
delta_gr <- function(x,y, ... ,dv=2){
 results<-orddom(x,y, ...)
 if (is.vector(x)) { x <- t(matrix(sort(x),1)) }
 if (is.vector(y)) { y <- t(matrix(sort(y),1)) }
 if (is.list(x)||is.data.frame(x)) { x <- as.matrix(x) }
 if (is.list(y)||is.data.frame(y)) { y <- as.matrix(y) }
 dv<-as.integer(dv)
if (ncol(results)==2){dv<-1}#which values in unpaired case
if((dv<1)||(dv>3)) {stop("the dv value of the delta_gr function must be between 1 and 3")}
results2<-as.matrix(as.integer(as.numeric(results[4:nrow(results),dv])*100000+.5)/100000)
rownames(results2)<-rownames(results)[4:nrow(results)]
if((results['H1 tails p/CI',dv]==1)&&(results['delta',dv]!=0)){#one-tailed
if(results['delta',dv]<0){#CI high
plotdata<-cbind(c(results2['delta',1],results2['delta',1],results2['CI high',1]))
}else{#CI low
plotdata<-cbind(c(results2['CI low',1],results2['delta',1],results2['delta',1]))
}}else{#two-tailed
plotdata<-cbind(c(results2['CI low',1],results2['delta',1],results2['CI high',1]))
if(dv==3) {plotdata<-cbind(c(results2['delta',1],results2['delta',1],results2['delta',1]))}
}
probn<-c(results2['N #Y<X',1],results2['N #Y=X',1],results2['N #Y>X',1])
nocomp<-sum(probn)
probn<-as.integer(probn/nocomp*10000+.5)/100
height<-250
#labels and titles for output
title<-"Cliff's delta "
label<-""
if(ncol(results)==4){
title<-paste(title," (",colnames(results)[dv],") ",sep="")
label<-paste("(",substr(colnames(results)[dv],1,1),")",sep="")}
#graphics output delta and ci
plot(cbind(plotdata,c(-15,-15,-15)),ylab="",xlab="",main=paste(title," & ",results['1-alpha',dv],"% Confidence Interval (",results['H1 tails p/CI',dv],"-tailed)",sep=""),type="b",xlim=c(-1,1),xaxp=c(-1,1,8),xaxs="i",ylim=c(-20,height),yaxt="n",lty=2,yaxs="i")
lines(x=c(plotdata[1,1],plotdata[1,1]),y=c(-17,-13),lty=1)
lines(x=c(plotdata[3,1],plotdata[3,1]),y=c(-17,-13),lty=1)
text(results['delta',dv],-15,pos=3,labels=paste("d",label,sep=""))
lines(x=c(0,0),y=c(-20,0),lty=3)
lines(x=c(-1,1),y=c(100,100),lty=1)
#output base
barplot(probn,beside=TRUE,space=c(-1.5,0,0),width=.6666666,ylim=c(0,119),yaxp=c(0,100,5),yaxs="i",col=c(grey(1),grey(0),grey(.5)),add=TRUE)
mtext(paste("% of all ordinal comparisons (",nocomp,")",sep=""),side=2,adj=0.125,line=2)
text(x=-(2/3),y=probn[1]+2,labels=paste(as.integer(probn[1]*100+.5)/100,"% (Y<X)\n",sep=""))
text(x=0,y=probn[2]+2,labels=paste(as.integer(probn[2]*100+.5)/100,"% (Y=X)\n",sep=""))
text(x=(2/3),y=probn[3]+2,labels=paste(as.integer(probn[3]*100+.5)/100,"% (Y>X)\n",sep=""))
#prepare legend
if (ncol(results)==2) {
 textlegende<-c(paste("X: ",results['var1_X',1]," (n=",length(x),")",sep=""))
 textlegende<-append(textlegende,paste("Y: ",results['var2_Y',1]," (n=",length(y),")",sep=""))
} else {
 textlegende<-c(paste("X: ",results['var1_X_pre',1]," (n=",length(x),")",sep=""))
 textlegende<-append(textlegende,paste("Y: ",results['var2_Y_post',1]," (n=",length(y),")",sep=""))
}
textlegende<-append(textlegende,paste("d",label," = ",results2['delta',1],sep=""))
textlegende<-append(textlegende,paste("s",label," = ",results2['s delta',1],sep=""))
textlegende<-append(textlegende,paste("z",label," = ",results2['z/t score',1],sep=""))
textlegende<-append(textlegende,paste("p",label," = ",results2['p',1],sep=""))
textlegende<-append(textlegende,paste("Cohen's d = ",results2["Cohen's d",1],sep=""))
lines(x=c(0,0),y=c(100,185),lty=1)
lines(x=c(-1,1),y=c(185,185),lty=1)
#data histograms
basepoint<-112
toppoint<-170
cohd<-as.numeric(results["Cohen's d",1])
sides<-sign(cohd)
cohd<-cohd*.10/2*sides
nx<- ((x-mean(x))/as.numeric((sqrt(var(x))))*.1)+.5-(cohd*sides)
ny<- ((y-mean(y))/as.numeric((sqrt(var(y))))*.1)+.5+(cohd*sides)
mx<-hist(nx,plot=FALSE)
mx$counts<-mx$counts/sum(mx$counts)
my<-hist(ny,plot=FALSE)
my$counts<-my$counts/sum(my$counts)
highest<-max(c(mx$counts,my$counts))
mx$counts<-mx$counts/highest*(toppoint-basepoint)
mx$counts<-mx$counts+basepoint
my$counts<-my$counts/highest*(toppoint-basepoint)
my$counts<-my$counts+basepoint
clip(.05,.95,basepoint,170)
plot(mx,add=TRUE,xlim=c(.10,90),xpd=FALSE)
plot(my,add=TRUE,xlim=c(.10,90),xpd=FALSE,col=grey(.5))
clip(-2,2,-2,height+20)
#overlap output
lines(x=c(.05,.95),y=c(basepoint,basepoint),lty=1)
first <- function(x) (((toppoint-basepoint)/dnorm(.5,.5,.1))+1)*dnorm(x,(.50-cohd),.10)+basepoint
second <- function(x) (((toppoint-basepoint)/dnorm(.5,.5,.1))+1)*dnorm(x,(.50+cohd),.10)+basepoint
plot(first,.10,.90,add=TRUE)
text(x=(.50-(cohd*sides)),y=toppoint,pos=3,labels="X")
plot(second,.10,.90,add=TRUE)
text(x=(.50+(sides*cohd)),y=toppoint,pos=3,labels="Y")
xt <- c(.05,seq(.05,.95,0.01),.95) 
yt <- c(basepoint,((((toppoint-basepoint)/dnorm(.5,.5,.1))+1)*dnorm(seq(.05,.5,0.01),(.50+cohd),.10)+basepoint),((((toppoint-basepoint)/dnorm(.5,.5,.1))+1)*dnorm(seq(.51,.95,0.01),(.50-cohd),.10)+basepoint),basepoint) 
polygon(xt,yt,col=rgb(.40,.40,.40,alpha=.33))
lines(x=c((.50-cohd),(.50-cohd)),y=c(basepoint,toppoint),lty=2)
lines(x=c((.50+cohd),(.50+cohd)),y=c(basepoint,toppoint),lty=2)
#output legend
legend(x=-1,y=188,cex=.9,legend=textlegende,bty="n",fill=c(grey(1),grey(.5),grey(1),grey(1),grey(1),grey(1)),border=c(grey(0),grey(0),grey(1),grey(1),grey(1),grey(1),grey(1)))
#interpretation of results
intprt<-""
if (ncol(results)==2){#unpaired case-----------------------------
intprt<-paste("There is a ",probn[3],"% chance that a case randomly chosen from group Y (",results['var2_Y',1],")\nis has a higher score than a randomly chosen subject or case from group X\n",sep="")
intprt<-paste(intprt,"(",results['var1_X',1],", as compared to a ",probn[1],"% probability for the reverse, resulting in a\nsuccess rate difference (i.e. ",sep="")
intprt<-paste(intprt,"an estimated delta value of) ",results2['delta',1],".\n",sep="")
} else {#paired case---------------------------------------------
if(dv==1){#within
intprt<-paste("Pairwise comparison yield higher scores in/at Y (",results['var2_Y_post',dv],")\nin ",probn[3],"% of all cases, as compared to a ",probn[1],"% probability ",sep="")
intprt<-paste(intprt,"to find higher scores in/at\nX (",results['var1_X_pre',1],"), resulting in an estimated delta (w) value of ",results2['delta',1],".\n",sep="")
}
if(dv==2){#between
intprt<-paste(intprt,"Excluding all ",length(x)," pairwise comparisons, higher scores in/at Y (",results['var2_Y_post',dv],")\nare found in ",probn[3],"% of cases when compared to all scores in/at\n",sep="")
intprt<-paste(intprt,"X (",results['var1_X_pre',1],"), resulting in an estimated 'group-shift' or\ndelta (b) value of ",results2['delta',1],".\n",sep="")
}
if(dv==3){#combined
intprt<-paste("The proportion of scores that are higher in/at Y (",results['var2_Y_post',1],")\nis ",probn[3],"%, as compared to a ",probn[1],"% probability ",sep="")
intprt<-paste(intprt,"to find higher scores in/at\nX (",results['var1_X_pre',1],"). Adding the effect of pairwise ordinal comparison and the\n",sep="")
intprt<-paste(intprt,"total ordinal group 'shift' yields an estimated combined\ndelta value of ",results2['delta',1],".\n",sep="")
}}
#output interpretation
intprt<-paste(intprt,"There is a probability of ",results['1-alpha',1],"% that the 'true' delta value ",sep="")
if((results['H1 tails p/CI',1]==1)&&(results['delta',1]!=0)){#one-tailed
if(results['delta',1]<0){#CI high
intprt<-paste(intprt,"is lower than\nCI(upper)= ",results2['CI high',1],".",sep="")
}else{#CI low
intprt<-paste(intprt,"is higher than\nCI(lower)= ",results2['CI low',1],".",sep="")
}}else{#two-tailed
if(results['H1 tails p/CI',1]==1){
intprt<-paste(intprt,"is either higher than\nCI(lower)= ",results2['CI low',1]," or lower than CI(higher)= ",results2['CI high',1],".",sep="")} else{
intprt<-paste(intprt,"can be found between\nCI(lower)= ",results2['CI low',1]," and CI(higher)= ",results2['CI high',1],".",sep="")}}
text(x=-.98,y=height-2,adj=c(0,1),cex=.8,labels=intprt)
}

