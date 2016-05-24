rsbplot<-function(data,plot.pval="TRUE",ylim.scan=2,pch=16,main=NA){
  tmp.plot=0
  if(colnames(data)[4]=="Pvalue (bilateral)"){yleg=expression("-" * log[10] * "[" ~ "1-2|" * Phi[scriptstyle(italic(rSB))] * "-0.5|" ~ "]");tmp.plot=1}
  if(colnames(data)[4]=="Pvalue (unilateral)"){yleg=expression("-" * log[10] * "[" ~ "1-" * Phi[scriptstyle(italic(rSB))] * "" ~ "]");tmp.plot=1}
  if(tmp.plot==0){
   warning("Unrecognized column name for Pvalue: plot.pval has been turned off")
   plot.pval="FALSE"
    }
    
   if(is.na(main)){main=colnames(data)[3]}

#   if(plot.pval){layout(matrix(1:2,2,1))}
   tmp_chr=unique(data[,1]) ; col_chr=1:length(tmp_chr) ; names(col_chr)=tmp_chr ;pos_chr=rep(0,length(tmp_chr))
   tmp_nmrk=table(data[,1]) ; pos_mrk=cumsum(tmp_nmrk)
   pos_chr[1]<-floor(pos_mrk[1]/2)
   if(length(tmp_chr)>1){for (i in 2:length(tmp_chr)){pos_chr[i]=pos_mrk[i-1] + floor((tmp_nmrk[i]/2))}}
   
   plot(data[,3],pch=pch,las=1,col=col_chr[as.character(data[,1])],xaxt="n",xlab="Chromosome",main=main,ylab="rSB")
   abline(h=ylim.scan,lty=2) ; abline(h=-1*ylim.scan,lty=2) 
   axis(1,at=pos_chr,labels=tmp_chr,las=1)
   
   if(plot.pval){
     plot(data[,4],pch=pch,las=1,col=col_chr[as.character(data[,1])],xaxt="n",xlab="Chromosome",main="Pvalue",
           ylab=yleg) 
      abline(h=ylim.scan,lty=2) ; abline(h=-1*ylim.scan,lty=2) 
     axis(1,at=pos_chr,labels=tmp_chr,las=1)
     abline(h=ylim.scan,lty=2) ; abline(h=-1*ylim.scan,lty=2) 
                }
#layout(matrix(1,1,1))
  }
