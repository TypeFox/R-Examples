`retrouve.all` <-
function(assoc,res.structure,all.marq){

     nbre.marq	=	length(all.marq)
     Hap	=	res.structure$Hap
     num.all	=	Hap[assoc,]
     all	=	rep(NA,nbre.marq)

     for(i in 1:nbre.marq){

         all[i]=all.marq[[i]][num.all[i]]

     }

     res    =	paste(all,collapse="")

     res

}

