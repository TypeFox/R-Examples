calc_ehhs<-function(haplohh,mrk,limhaplo=2,limehhs=0.05,plotehhs=TRUE,main_leg="EHHS plot"){
  if(!(is.haplohh(haplohh))){stop("Data oject is not of valid haplohh object... (see data2haplohh() function)")} 
  if(mrk<1 | mrk>haplohh@nsnp){stop(paste("Focal snp index must be between",1,"and",haplohh@nsnp))}
  if(limhaplo<2){stop("limhaplo must be >1")}
  if(limehhs<0 | limehhs>1){stop("limehhs must be between 0 and 1")}
  
  nhaplo_eval<-rep(0,haplohh@nsnp) ; ehhs<-rep(0,haplohh@nsnp) ; ies<-0
  res.ehhs<-.C("r_ehhs", 
                  Rdata = as.integer(haplohh@haplo),
                  number_SNPs  = as.integer(haplohh@nsnp),
                  number_chromosomes = as.integer(haplohh@nhap),
                  focal_SNP = as.integer(mrk),
                  map = as.double(haplohh@position),
                  number_haplotypes = as.integer(nhaplo_eval),
                  EHHS = as.double(ehhs),
                  IES = as.double(ies),
                  min_number_haplotypes = as.integer(limhaplo),
                  min_EHH = as.double(limehhs)
                  )

  ehhs=res.ehhs$EHHS ; nhaplo_eval=res.ehhs$number_haplotypes 
  names(ehhs)=names(nhaplo_eval)=haplohh@snp.name
  
  
 if(plotehhs){
   sel_reg<-(nhaplo_eval>0)
   if(sum(sel_reg)>0){
     plot(haplohh@position[sel_reg],ehhs[sel_reg],col=c("red"),lty=1,
             type="l",main=main_leg,bty="n",xlab="Position",ylab="EHHS")
     abline(v=haplohh@position[mrk],lty=2)
   }
   }
   
  return(list(ehhs=ehhs,nhaplo_eval=nhaplo_eval,ies=res.ehhs$IES))
}

