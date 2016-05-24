`recode.hap` <-
function(hap,all.marq){


    nb.marq	=	length(hap[1,])
    nb.ind	=	length(hap[,1])

    hap.recode=	matrix(NA,ncol=nb.marq,nrow=nb.ind)

    for(i in 1:nb.marq){

         for(k in 1:length(all.marq[[i]])){

             hap.recode[hap[,i]==all.marq[[i]][k],i]=k
         }
     }


hap.recode


}

