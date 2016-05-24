`corresp` <-
function(hap.o,res.structure){

hap.obs	=	as.matrix(hap.o)
Bin.type=	res.structure$Bin.type
Hap	=	res.structure$Hap
nbre.type=	length(Bin.type[,1])
nbre.hap=	length(Hap[,1])
nbre.ind=	length(hap.obs[,1])
nbre.marq=	length(Bin.type[1,])
corresp		=	matrix(NA,nrow=nbre.ind,ncol=4)
corresp[,1]	=	1:nbre.ind

#création des haplotypes possibles
nbre.hap.poss	=	rep(0,nbre.type)
hap.poss	=	list()

   for(ik in (2:nbre.type)){
       iik=nbre.type-(ik-1)
       Bin=Bin.type[ik,]
       One=(1:nbre.marq)[Bin==1]
       hap.colle=rep(NA,nbre.hap)

       for (ih in 1:nbre.hap){
            hap.colle[ih]=paste(Hap[ih,One],collapse="")
       }
       hap.poss[[iik]]=unique(hap.colle)
       nbre.hap.poss[iik]=length(hap.poss[[iik]])
    }


exposant	=	(nbre.marq-1):0

   for(ind in 1:nbre.ind){

       ABS=(1:nbre.marq)[is.na(hap.obs[ind,])]
       PRE=(1:nbre.marq)[!is.na(hap.obs[ind,])]
       bi.type=rep(0,nbre.marq)

       if(length(ABS) !=0) bi.type[ABS]=1
       id.type=sum(bi.type*2^exposant)+1
       id.hap=1

       if(length(PRE) !=0) {
          hap.obs.colle=paste(hap.obs[ind,PRE],collapse="")
          id.hap=(1:length(hap.poss[[id.type]]))[hap.poss[[id.type]]==hap.obs.colle]
       }

       corresp[ind,2]	=	id.type
       corresp[ind,3]	=	id.hap
    }


corresp[,4]	=	corresp[,3]
which		=	(1:nbre.ind)[corresp[,2]>1]

   for(iw in which)  {
      corresp[iw,4]=corresp[iw,4]+sum(nbre.hap.poss[1:(corresp[iw,2]-1)])
   }


assoc	=	unique(corresp[corresp[,2]==1,3])

res	=	list("corresp"=corresp,"assoc"=assoc)


res

}

