plot.NAM = function(x,...,alpha=0.05,colA=2,colB=4,find=NULL,FDR=NULL,gtz=FALSE,phys=NULL){
  anyNA = function(x) any(is.na(x))
  if(!is.null(FDR)){if(FDR>=1|FDR<0)stop("FRD must be between 0 and 1")}
  gwas=x
  chr=as.numeric(summary(factor(as.numeric(gwas$MAP[,1]))))
  
  if(gwas$Method=="P3D"){
    
    FGWASplot=function(Fgwas,chr,AA,BB,...){
      col=c();for(i in 1:length(chr)){if((i%%2)==0){b=AA}else{b=BB};a=rep(b,chr[i]);col=c(col,a)}
      W = Fgwas$PolyTest$wald
      if(is.null(phys)){ Xaxis = 1:length(W) }else{ Xaxis = cumsum(phys) }
      plot(Xaxis,W,col=col,xlab="Chromosome",ylab="Wald Statistics",xaxt = "n",...)
      return(W)}
    
    # Plot
    pv=FGWASplot(gwas,chr=chr,AA=colA,BB=colB,...)
    
    # QTL
    if(!is.null(find)){Loc=identify(n=find,x=1:length(pv),y=pv,labels=gwas$SNPs);for(i in Loc) cat(gwas$SNPs[i],'\n')}

    }else{
    
    RGWASplot=function(Rgwas,chr,AA,BB,...){
      # Colors
      col=c();for(i in 1:length(chr)){if((i%%2)==0){b=AA}else{b=BB};a=rep(b,chr[i]); col=c(col,a)}
      # Statistics
      S = Rgwas$PolyTest$pval
      if(is.null(phys)){ Xaxis = 1:length(S) }else{ Xaxis = cumsum(phys) }
      plot(Xaxis,S,col=col,xlab="Chromosome",ylab="-log(p-value)",xaxt = "n",...)
      return(S)
    }
    
    # Plot
    if (is.null(alpha)){
      
      pv=RGWASplot(gwas,chr=chr,AA=colA,BB=colB,...)
  
      } else {
      
      pv=RGWASplot(gwas,chr=chr,AA=colA,BB=colB,...)
    
      if(is.null(FDR)){
        A = 1-alpha
        LRmax = qchisq(A,0.5)
        lim = -log(dchisq(LRmax, 0.5),base = 10)
        abline(h=lim,col=1,lty=2) 
      }else{
        
        # Computing where chromosomes start and end
        NumChr = length(chr)
        Ch0 = cumsum(c(1,chr[-NumChr]))-.5
        Ch1 = cumsum(chr)+.5
       
        # Multiple test correction
        if(gtz==T){
          
          MT = tapply(
            X = gwas$PolyTest$pval,
            INDEX = gwas$MAP[,1],
            FUN = function(x) sum(x>0)
            )
          MT[MT==0]=1
          
        }else{MT=chr}
        
        
        # Thresholds
        for(i in 1:(NumChr)){
          A = 1-alpha/(MT[i]*(1-FDR))
          LRmax = qchisq(A,0.5)
          lim = -log(dchisq(LRmax, 0.5),base = 10)
          lines(x = c(Ch0[i],Ch1[i]),y = c(lim,lim))
        }
        
      }
      }
     
    # QTL
    if(!is.null(find)){Loc=identify(n=find,x=1:length(pv),y=pv,labels=gwas$SNPs);for(i in Loc) cat(gwas$SNPs[i],'\n')}
    
    
    }
  
  # Adding Chromosome in X axis
  medians=rep(NA,length(chr))
  
  if(is.null(phys)){
    for(i in 1:length(chr)) medians[i] = median(which(gwas$MAP[,1]==i))
    axis(1, at=round(medians), labels=1:length(medians))
  }else{
    Xaxis = cumsum(phys)
    #for(i in 1:length(chr)) medians[i] = median(Xaxis[gwas$MAP[,1]==i])
    for(i in 1:length(chr)) medians[i] = mean(range(Xaxis[gwas$MAP[,1]==i]))
    axis(1, at=round(medians), labels=1:length(medians))
  }
  
  
}