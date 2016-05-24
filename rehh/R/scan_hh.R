scan_hh<-function(haplohh,limhaplo=2,limehh=0.05,limehhs=0.05){
  if(!(is.haplohh(haplohh))){stop("Data oject is not of valid haplohh object... (see data2haplohh() function)")}
  if(limhaplo<2){stop("limhaplo must be >1")}
  if(limehh<0 | limehh>1){stop("limehh must be between 0 and 1")}
  if(limehhs<0 | limehhs>1){stop("limehhs must be between 0 and 1")}
    IHH <- matrix(0.0,nrow = haplohh@nsnp,ncol = 2)
    IES <- vector(mode = "numeric",length = haplohh@nsnp)
    res_scan<-.C("r_scan_hh", 
                    Rdata = as.integer(haplohh@haplo),
                    number_SNPs  = as.integer(haplohh@nsnp),
                    number_chromosomes = as.integer(haplohh@nhap),
                    map = as.double(haplohh@position),
                    IHH = as.double(IHH),
                    IES = as.double(IES),
                    min_number_haplotypes = as.integer(limhaplo),
                    min_EHH = as.double(limehh),
                    min_EHHS = as.double(limehhs)
                    )
    tmp_n1=colSums(haplohh@haplo==1) ; tmp_freq=tmp_n1/(tmp_n1 + colSums(haplohh@haplo==2))
    RES_ALL=cbind(rep(haplohh@chr.name,haplohh@nsnp),haplohh@position,tmp_freq,matrix(res_scan$IHH,haplohh@nsnp,2),res_scan$IES)
    rownames(RES_ALL)=haplohh@snp.name ; colnames(RES_ALL)=c("CHR","POSITION","FREQ_a","IHHa","IHHd","IES")
  return(RES_ALL)

}
