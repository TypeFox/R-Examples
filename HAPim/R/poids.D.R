`poids.D` <-
function(dist.test,pos.QTL,res.structure){

    Bin.type	=	res.structure$Bin.type
    nbre.type	=	length(Bin.type[,1])
    nbre.marq	=	length(dist.test)-1
    temp	=	rep(NA,(nbre.marq+1))
    poids.D	=	list()
    poids.D[[1]]=	1

   for(ik in 2:(nbre.type)){

      Bin                =    Bin.type[ik,]
      temp[1:(pos.QTL-1)]=    Bin[1:(pos.QTL-1)]
      temp[pos.QTL]      =    1
      temp[(pos.QTL+1):(nbre.marq+1)]=Bin[pos.QTL:nbre.marq]
      dist.temp          =    dist.test[temp==1]
      Ord                =    length(dist.temp)
      d1        =      dist.temp[2:Ord]
      d2        =      dist.temp[1:(Ord-1)]
      dist.temp =      d1-d2

      poids.D[[ik]]=1-haldanem1(dist.temp)

   }

poids.D


}

