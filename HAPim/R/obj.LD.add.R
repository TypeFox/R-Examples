`obj.LD.add` <-
function(param,don)
{
 perf            = don[[1]]  
 CD              = don[[2]]
 assoc           = don[[3]]
 cor.pere        = don[[4]]
 cor.mere        = don[[5]]
 pi.hap          = don[[6]]
 res.structure   = don[[7]]
 poids.D         = don[[8]]
 nind            = length(perf)

 if(param[1] >0 & param[3]>0 & param[4]>0  & param[4]<1 ) {
  
 timeT=param[3]
 piQ.t0=param[4]
  
 esp.freq.hap=esp.freq.hap(assoc,piQ.t0,timeT,pi.hap,res.structure,poids.D)
 DL.pere=proba.DL(piQ.t0,esp.freq.hap,res.structure,pi.hap,cor.pere)
 DL.mere=proba.DL(piQ.t0,esp.freq.hap,res.structure,pi.hap,cor.mere)

 DL.d=proba.DL.diplotype(DL.pere,DL.mere)
 vraisemb=-vrais.LD.add(param[5],param[2],param[1],CD,perf,DL.d) 

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

