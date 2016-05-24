`proba.DL` <-
function(piQ.t0,esp.freq.hap,res.structure,pi.hap,res.corresp){

    corresp   =    res.corresp$corresp
    Bin.type  = res.structure$Bin.type
    Index     = res.structure$Index
    nbre.type = length(Bin.type[,1])
    nbre.ind  = length(corresp[,1])
    DL        = list()
    DL[[1]]   = esp.freq.hap/pi.hap[[1]]
    DL[[1]][DL[[1]]>1]=1
    DL[[nbre.type]]=piQ.t0

    for(ik in 2:(nbre.type-1)){
       Ind=Index[[ik]]
       iik=nbre.type-(ik-1)
       nb.h=length(Ind[,1])
       DL[[iik]]=rep(NA,nb.h)
       temp=esp.freq.hap/pi.hap[[iik]]

           for(ih in 1:nb.h){
               DL[[iik]][ih]=sum(temp[Ind[ih,]])
            }
     DL[[iik]][DL[[iik]]>1]=1
     }

    DL.vect=NULL

   for(ik in 1:nbre.type){
   DL.vect=c(DL.vect,DL[[ik]])
   }  

proba.DL=DL.vect[corresp[,4]]

proba.DL

}

