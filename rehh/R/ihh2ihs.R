ihh2ihs<-function(res_ihh,freqbin=0.025,minmaf=0.05){
  res_ihh=res_ihh[res_ihh[,3]>minmaf & res_ihh[,3]<(1-minmaf) , ]
  res_ihs=cbind(res_ihh[,1:2],rep(0,nrow(res_ihh)))
  colnames(res_ihs)[3]="iHS" ; rownames(res_ihs)=rownames(res_ihh)
  if(freqbin>0){
    freq_class=seq(minmaf,1-minmaf,freqbin)
    summary_class=matrix(0,length(freq_class)-1,4) ; colnames(summary_class)=c("Freq Class","Size","Mean iHH","SD iHH")
    ihs=log(res_ihh[,4]/res_ihh[,5]) ; ihs[ihs=="Inf" | ihs=="-Inf"]=NA 
     for(c in 1:(length(freq_class)-1)){
       lim_inf=freq_class[c] ; lim_sup=freq_class[c+1]
       mrk_sel=(res_ihh[,3]>=lim_inf & res_ihh[,3]<lim_sup)
       tmp_ihs=ihs[mrk_sel] ; tmp_nmrk=sum(mrk_sel)
       if(tmp_nmrk<10){
        warning(paste("Size of Allele Frequency Class: ",lim_inf,"-",lim_sup," <10: You should probably increase freqbin\n",sep=""))
           }
       tmp_mean=mean(tmp_ihs,na.rm=T) ; tmp_sd=sd(tmp_ihs,na.rm=T)
       summary_class[c,1]=paste(lim_inf,"-",lim_sup) ; summary_class[c,2]=tmp_nmrk 
       summary_class[c,3]=tmp_mean ;  summary_class[c,4]=tmp_sd
       ihs[mrk_sel]=(ihs[mrk_sel]-tmp_mean)/tmp_sd
            }
  }else{
    freq_class=unique(res_ihh[,3]) ; freq_class=freq_class[order(freq_class)]
    summary_class=matrix(0,length(freq_class),4) ; colnames(summary_class)=c("Freq Class","Size","Mean iHH","SD iHH")
    ihs=log(res_ihh[,4]/res_ihh[,5]) ; tmp_std=matrix(NA, length(ihs), 3) #NUM,MEAN,SD
    for(f in 1:length(freq_class)){
      mrk_sel=(res_ihh[,3]==freq_class[f]) ; tmp_ihs=ihs[mrk_sel]
      tmp_mean=mean(tmp_ihs,na.rm=T) ; tmp_sd=sd(tmp_ihs,na.rm=T)
            summary_class[f, 1] = freq_class[f] ; summary_class[f, 2] = sum(mrk_sel) 
            summary_class[f, 3] = tmp_mean ; summary_class[f, 4] = tmp_sd  
      ihs[mrk_sel]=(ihs[mrk_sel]-tmp_mean)/tmp_sd
 }
  }
  res_ihs[,3]=ihs ; tmp_pval=-1*log10(1-2*abs(pnorm(ihs)-0.5))
  tmp_pval2=tmp_pval ; tmp_pval2[tmp_pval2=="Inf"]=NA 
  tmp_pval[tmp_pval=="Inf"]=max(tmp_pval2,na.rm=TRUE) + 1 
  res_ihs=cbind(res_ihs,tmp_pval) ; colnames(res_ihs)[4]="Pvalue"

  return(list(res.ihs=res_ihs,summary.class=summary_class))
}
