`degenset` <-
function(data){

 gooditem <- colnames(data)
 gooditem <- gooditem[apply(data[,gooditem],2,function(XXX){temp <- unique(XXX)
                                       temp <- temp[! is.na(temp)]
                                       length(temp) > 1})]
 
 repeat{
  u.data     <- unique.data.counts(data[,gooditem])
  gooditem   <- colnames(u.data$x.i)
  goodchange <- length(gooditem)
  gooditem   <- gooditem[apply(u.data$x.i[,gooditem],2,function(XXX){temp <- unique(XXX)
                                       temp <- temp[! is.na(temp)]
                                       length(temp) > 1})]
  if(length(gooditem) == goodchange) break else goodchange <- length(gooditem)
 }

 baditems <- colnames(data)[! colnames(data) %in% gooditem] 
 if(length(baditems) > 0) warning(paste("Removed items with perfect scores: ",baditems,"\n",sep=""))
 
 list(u.data,baditems)
}

