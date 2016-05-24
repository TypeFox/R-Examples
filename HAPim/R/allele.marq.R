`allele.marq` <-
function(hap){

nb.marq		=	length(hap[1,])
all.marq	=	list()

     for(i in 1:nb.marq){

          all.marq[[i]]=unique(na.omit(hap[,i]))
      }


all.marq



}

