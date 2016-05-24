GWplot<-function(data,pval,pos,chr,chr.plot=c(1:22,"X"),title.text="",ylim=Inf,outfile,cutoff1=5e-8,cutoff2=4e-7){
        data<-data[!is.na(data[,pval])&!is.na(data[,pos])& data[,chr]%in%chr.plot,]
        ##change chrosome that are not 1:22 or X to Other
        chrs<-c(1:length(unique(data[,chr])))
        nchr<-length(chrs)
        bitmap(outfile,width=12,height=8)
        par(mar=c(5,5,4,2))
        if (any(data[,pval]==0)){
            warning("0 p-values changed to 5E-324, see pvalzero....csv")
            write.table(data[data[,pval]==0,],paste("pvalzero.csv"),sep=",")
            data[data[,pval]==0,pval]<-5E-324    
        }
        value<-names(table(data[,chr]))[table(data[,chr])!=0]
        #print(value)
        if (length(value)>22){
           for(i in 23:nchr){
               #data[data[,chr]==value[i],chr]=i
               chrom <- rep(0,nrow(data))
               chrom[data[,chr]==value[i]] <- i
               chrom[chrom!=i] <- as.numeric(as.character(data[data[,chr]!=value[i],chr]))
               data[,chr] <- chrom
           }
        }
        data[,chr]<-as.numeric(data[,chr])
        phy.max<-tapply(data[,pos], data[,chr],max,na.rm=T)
        cumlen=0
        for(i in 1:nchr){
           data[data[,chr]==i,"loc"]<-data[data[,chr]==i,pos]+cumlen
           cumlen<-cumlen+phy.max[i]
        }
        phy.med<-tapply(data[,"loc"], data[,chr],median,na.rm=T)
        #print(phy.med)    
        data$mlgpval<--log(data[,pval], base=10)
        #print(max(data$mlgpval,na.rm=T))
        #print(min(max(data$mlgpval,na.rm=T)+2,Inf))
 	 if (ylim==Inf) ylim=min(max(data$mlgpval,na.rm=T)+2,Inf)
        plot(data[,"loc"],data[,"mlgpval"],type="n",xaxt="n",frame.plot=F,
             main=title.text, cex=0.5,cex.axis=1.5,cex.lab=1.5,xlab="",ylab="",
             xlim=c(0,max(data[,"loc"],na.rm=T)),ylim=c(0,ylim))
        if (nchr==22) label=1:22 else label=c(1:22,value[23:nchr])
	 axis(side=1, at=phy.med, labels=label,cex.lab=0.8)
        for(i in 1:nchr){
           if(i %in% seq(2,30,2)) col="grey" else col="darkgrey"
           points(data[data[,chr]==i&data$mlgpval<(-log10(cutoff2)),"loc"],
           data[data[,chr]==i&data$mlgpval<(-log10(cutoff2)),"mlgpval"],col=col,pch=20,cex=0.6)
           points(data[data[,chr]==i&data$mlgpval<(-log10(cutoff1))&data$mlgpval>=(-log10(cutoff2)),"loc"],
           data[data[,chr]==i&data$mlgpval<(-log10(cutoff1))&data$mlgpval>=(-log10(cutoff2)),"mlgpval"],col="blue",pch=20,cex=0.6)
           points(data[data[,chr]==i&data$mlgpval>=(-log10(cutoff1)),"loc"],
           data[data[,chr]==i&data$mlgpval>=(-log10(cutoff1)),"mlgpval"],col="red",pch=20,cex=0.6)
           #abline(4,0,col="grey")
           abline(-log10(cutoff1),0,col="grey",lty=2)
        }
        dev.off()
}

