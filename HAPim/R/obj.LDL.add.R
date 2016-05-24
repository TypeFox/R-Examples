`obj.LDL.add` <-
function(param,don)
{
 perf            = don[[1]]  
 CD              = don[[2]]
 PLA             = don[[3]]
 desc.pere       = don[[4]]
 moyenne.pere    = don[[5]]
 assoc           = don[[6]]
 cor.chrm1.pere  = don[[7]]
 cor.chrm2.pere  = don[[8]]
 cor.mere        = don[[9]]
 pi.hap          = don[[10]]
 res.structure   = don[[11]]
 poids.D         = don[[12]]
 nind            = length(perf)

 if(param[1] >0 & param[3]>0 & param[4]>0 & param[4]<1) {

# on commence par les pères au temps t

  timeT  = param[3]
  piQ.t0 = param[4]
  
  esp.freq.hap=esp.freq.hap(assoc,piQ.t0,timeT,pi.hap,res.structure,poids.D)
  DL.chrom1=proba.DL(piQ.t0,esp.freq.hap,res.structure,pi.hap,cor.chrm1.pere)
  DL.chrom2=proba.DL(piQ.t0,esp.freq.hap,res.structure,pi.hap,cor.chrm2.pere)


# fin des pères


# début des mères au temps t+1

  esp.freq.hap=esp.freq.hap(assoc,piQ.t0,(timeT+1),pi.hap,res.structure,poids.D)
  DL.m=proba.DL(piQ.t0,esp.freq.hap,res.structure,pi.hap,cor.mere)

  vraisemb = -vrais.LDL.add(moyenne.pere,param[2],param[1],CD,perf,PLA,DL.m,DL.chrom1,DL.chrom2,desc.pere,param[5]) 

  }




  if(param[1]<0){ 
    param[1]=0
    vraisemb =750*nind
  }

  if(param[3]<0){
    param[3]=0
    vraisemb =750*nind
  }

  if(param[4]<0){
    param[4]=0
    vraisemb =750*nind
  }

  if(param[4]>1){
    param[4]=1
    vraisemb =750*nind
  }

 vraisemb  

         
}

