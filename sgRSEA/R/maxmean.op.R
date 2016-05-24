maxmean.op <-
function(zj){
     m=length(zj)
     sc.po =  sum( ifelse(zj>=0, zj, 0)) /m
     sc.ne =  sum( ifelse(zj<0, zj, 0) ) /m
     sc = c(sc.po, sc.ne)
     return(    tn=sc[which.max( abs(sc) ) ]  )
    }
