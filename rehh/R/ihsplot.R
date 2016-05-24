ihsplot<-function(data,plot.pval="TRUE",ylim.scan=2,pch=16,main="iHS"){
#   if(plot.pval){layout(matrix(1:2,2,1))}
   tmp_chr=unique(data[,1]) ; col_chr=1:length(tmp_chr) ; names(col_chr)=tmp_chr ;pos_chr=rep(0,length(tmp_chr))
   tmp_nmrk=table(data[,1]) ; pos_mrk=cumsum(tmp_nmrk)
   pos_chr[1]<-floor(pos_mrk[1]/2)
   if(length(tmp_chr)>1){for (i in 2:length(tmp_chr)){pos_chr[i]=pos_mrk[i-1] + floor((tmp_nmrk[i]/2))}}
   
   plot(data[,3],pch=pch,las=1,col=col_chr[as.character(data[,1])],xaxt="n",xlab="Chromosome",ylab="iHS",main=main)
   abline(h=ylim.scan,lty=2) ; abline(h=-1*ylim.scan,lty=2) 
   axis(1,at=pos_chr,labels=tmp_chr,las=1)
   
   if(plot.pval){
     plot(data[,4],pch=pch,las=1,col=col_chr[as.character(data[,1])],xaxt="n",xlab="Chromosome",main="Pvalue",
           ylab=expression("-" * log[10] * "[" ~ "1-2|" * Phi[scriptstyle(italic(iHS))] * "-0.5|" ~ "]")) 
      abline(h=ylim.scan,lty=2) ; abline(h=-1*ylim.scan,lty=2) 
     axis(1,at=pos_chr,labels=tmp_chr,las=1)
     abline(h=ylim.scan,lty=2) ; abline(h=-1*ylim.scan,lty=2) 
                }
#layout(matrix(1,1,1))
  }
