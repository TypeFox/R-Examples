#Plot method
setMethod("plot",
    signature(x = "MixMAP"),
   function(x,col.genes=c("black","gray"),col.detected=c("blue","violet"),col.text="black",title="MixMAP Manhattan Plot",display.text=TRUE){
############################
#defining errors
############################
#is input object of class MixMAP
if (class(x)!="MixMAP") stop(gettextf('input is not of class MixMAP'))

############################
#warnings
############################
if (length(col.genes)!=2) {(warning('\"col.genes\" should be a character vector of length 2 containing the names of colors  for alternating chromosomes'))}

#remove non-numeric chromosomes
manhat<-x@output[!grepl("[A-z]",x@output$chr),]

#sort the data by chromosome and then base pair
manhat.ord<-manhat[order(as.numeric(manhat$chr),manhat$coord),]
manhat.ord<-manhat.ord[!is.na(manhat.ord$coord),]


##Finding the maximum position for each chromosome
max.pos<-NULL
for (i in 1:21){
max.pos[i]<-max(manhat.ord$coord[manhat.ord$chr==i])}
max.pos1<-c(0,max.pos)
max.pos2<-NULL
for (i in 1:22){max.pos2[i]<-sum(max.pos1[1:i])}

#Add spacing between chromosomes
max.pos2<-max.pos2+c(0:21)*100000000

#defining the postitions of each gene in the plot
manhat.ord$pos<-manhat.ord$coord+max.pos2[as.numeric(manhat.ord$chr)]
manhat.ord$postEst[manhat.ord$postEst>0]<-0
manhat.ord$postEst[manhat.ord$postEst<0]=abs(manhat.ord$postEst[manhat.ord$postEst<0])

#defining the coloring for the MixManhattan plot
manhat.ord$col[as.numeric(manhat.ord$chr)%%2==0]<-col.genes[1]
manhat.ord$col[as.numeric(manhat.ord$chr)%%2==1]<-col.genes[2]

text.pos<-rep(NA,22)
for (i in 1:22){text.pos[i]<-mean(manhat.ord$pos[manhat.ord$chr==i])}

#plot the data
plot(manhat.ord$pos/1000000,manhat.ord$postEst,pch=20,col=manhat.ord$col,xlab="Chromosome",ylab="Absolute Value of Empirical Bayes Estimate",axes=F,main=title,ylim=c(0,max(manhat.ord$postEst)+0.5))
axis(2)
abline(h=0)

#which genes are detected in MixMAP only? Which genes are detected in both MixMAp and single SNP analysis
MixMAPGenes <- as.character(x@detected.genes[,1])
MixMAPandSNPGenes<-as.character(x@detected.genes[x@detected.genes$pval.min<(0.05/1000000),1])
MixMAPonlyGenes<-setdiff(MixMAPGenes,MixMAPandSNPGenes)

#Add legend
legend("topright",c("MixMAP positive; Single SNP negative","MixMAP positive; Single SNP positive"),border=col.detected,col=col.detected,pch=c(15,19),pt.cex=c(2,1))

#Add chromosome number
text(text.pos/1000000,-.05,seq(1,22,by=1),xpd=TRUE,cex=1)

#Plotting detected genes
#Were any genes detected?
if (x@num.genes.detected[1]>0){
#Plot the detected genes with a different symbol and color in MixMAP but not single SNP
points(manhat.ord$pos[manhat.ord[,1]%in%MixMAPonlyGenes]/1000000,manhat.ord$postEst[manhat.ord[,1]%in%MixMAPonlyGenes],pch=15,col=col.detected[1], bg=col.detected[1],cex=2)

#Were any genes detected in both SNP and MixMAP?
if (any(x@detected.genes$pval.min<0.05/1000000)){
#Plot MixMAP and Single SNP detected genes
points(manhat.ord$pos[manhat.ord[,1]%in%MixMAPandSNPGenes]/1000000,manhat.ord$postEst[manhat.ord[,1]%in%MixMAPandSNPGenes],pch=19,col=col.detected[2], bg=col.detected[2])
}

}

#Display gene names of detected genes
if (display.text==TRUE & x@num.genes.detected[1]>0 & length(MixMAPonlyGenes)>0){
	text(manhat.ord$pos[manhat.ord[,1]%in%MixMAPonlyGenes]/1000000,manhat.ord$postEst[manhat.ord[,1]%in%MixMAPonlyGenes],as.character(manhat.ord[manhat.ord[,1]%in%MixMAPonlyGenes,1]),col=col.text,offset=.6,pos=4, cex=.8)
}

if (display.text==TRUE & x@num.genes.detected[1]>0 & length(MixMAPandSNPGenes)>0){
text(manhat.ord$pos[manhat.ord[,1]%in%MixMAPandSNPGenes]/1000000,manhat.ord$postEst[manhat.ord[,1]%in%MixMAPandSNPGenes],as.character(manhat.ord[manhat.ord[,1]%in%MixMAPandSNPGenes,1]),col=col.text,offset=.6,pos=4, cex=.8)
}


})

