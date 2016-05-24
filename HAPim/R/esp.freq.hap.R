`esp.freq.hap` <-
function(hap.assoc,piQ.t0,timeT,pi.hap,res.structure,poids.D){

      Bin.type	=	res.structure$Bin.type
      Index	=	res.structure$Index
      Som.freq	=	res.structure$Som.freq
      nbre.marq	=	length(Bin.type[1,])
      nbre.type	=	length(Bin.type[,1])
      nbre.hap	=	length(pi.hap[[1]])

      Des.t0	=	matrix(NA,nrow=nbre.hap,ncol=nbre.type)
      Des.T	=	matrix(NA,nrow=nbre.hap,ncol=nbre.type)
      freq	=	matrix(1,nrow=nbre.hap,ncol=nbre.type)
      pi.Qh	=	rep(0,nbre.hap)
      pi.Qh[hap.assoc]=	piQ.t0

# frequence d'ordre 0, indépendance
     freq[,1] = rep(piQ.t0,nbre.hap)*pi.hap[[1]]

#fréquence d'ordre nbre.type, dépendance totale
     freq[,nbre.type] = pi.Qh

     for(ik in (2:(nbre.type-1))){

         Ind=Index[[ik]]
         temp=rep(0,nbre.hap)
         nbr.lig=nrow(Ind)

           for(nl in (1:nbr.lig)){

               which=Ind[nl,]
               temp[which]=sum(pi.Qh[which])

           }

         freq[,ik]=temp*pi.hap[[ik]]

     }

Des.t0	=	freq %*% Som.freq

#calcul des déséquilibres au temps timeT

     Des.T[,1]=Des.t0[,1]

     for(ik in  2:(nbre.type)){

        poids = poids.D[[ik]]
        poids = prod(poids^timeT)
        Des.T[,ik]=poids*Des.t0[,ik]

     }

#calcul de l'espérence des féquences des haplotypes et de Q au temps timeT
        Esp.freq.hap=rowSums(Des.T)
        Esp.freq.hap



}

