ies2rsb<-function(hh_pop1,hh_pop2,popname1=NA,popname2=NA,method="bilateral"){

  ies_1=hh_pop1[,6] ; ies_2=hh_pop2[,6]
  if(!(nrow(hh_pop1)==nrow(hh_pop2))){stop("hh_pop1 and hh_pop2 must have the same dimensions")}
  if(sum(hh_pop1[,2]==hh_pop2[,2])<nrow(hh_pop1)){stop("SNP position in hh_pop1 and hh_pop2 must be the same")}
  tmp_rsbnc=log(ies_1/ies_2) ; tmp_med=median(tmp_rsbnc,na.rm=T) ; tmp_sd=sd(tmp_rsbnc,na.rm=T)
  rsbcor=(tmp_rsbnc-tmp_med)/tmp_sd
  tmp_pval=rsbcor*0
  if(method=="bilateral"){tmp_pval=-1*log10(1-2*abs(pnorm(rsbcor)-0.5))}
  if(method=="unilateral"){tmp_pval=-1*log10(1-pnorm(rsbcor))}
  tmp_pval2=tmp_pval ; tmp_pval2[tmp_pval2=="Inf"]=NA  
  tmp_pval[tmp_pval=="Inf"]=max(tmp_pval2,na.rm=TRUE) + 1 
  res.rsb=cbind(hh_pop1[,1:2],rsbcor,tmp_pval)
  colnames(res.rsb)[3]=paste("rSB (",popname1," vs ",popname2,")",sep="")
  colnames(res.rsb)[4]=paste("Pvalue (",method,")",sep="")

  return(list(res.rsb=res.rsb))
}
