`pi.hap` <-
function(freq.marq,res.structure){

   Bin.type	=	res.structure$Bin.type
   Index	=	res.structure$Index
   Hap		=	res.structure$Hap
   nbre.type	=	length(Bin.type[,1])
   nbre.marq	=	length(freq.marq)
   nbre.hap	=	length(Hap[,1])
   pi.hap       =       list()
   pi.marq      =       matrix(NA,nrow=nbre.hap,ncol=nbre.marq)

   for(im in 1:nbre.marq){

       for(ia in 1:nbre.hap){

            which=(1:nbre.hap)[as.numeric(Hap[,im])==ia]
               pi.marq[which,im]=freq.marq[[im]][ia] 

        }

   }

   pi.hap[[nbre.type]]=rep(1,nbre.hap)

   for(it in 1:(nbre.type-1)){

       Bin=Bin.type[it,]
       Zero=(1:nbre.marq)[Bin==0]
       pi.hap[[it]]=apply(as.matrix(pi.marq[,Zero]),1,prod)

   }

pi.hap


}

