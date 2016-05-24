modChi2 <-
function(data,alp=0.5){
#Phase 1
        dat <- data
        p <- ncol(dat)-1      
        alpha <- alp; d <- 0.1 

        LevO <- LevCon(dat)
        Disc <- chiM(dat,alpha=alpha)$Disc.data
        LevN <- LevCon(Disc)
        e <- 0.01
     while((LevN-LevO)<e){
               if(alpha<=0.2) break
               alpha <- alpha-d
               Disc <- chiM(Disc,alpha=alpha)$Disc.data
               LevN <- LevCon(Disc) 
               LevO <- LevN
    }

#Phase 2
       options(digits=3)
       eps <- 0.0001
       sig <- array(alpha,p)
       d <- 10
    cutp=list()
    for(i in 1:p){
                LevO <- LevCon(Disc)
              while(TRUE){
                    val <- value(i,Disc,alpha=sig[i])
                    Disc <- val$disc
                    LevN <- LevCon(Disc)#inconRate for data
                    if(LevN<LevO|| sig[i]<=eps) break 
                    LevO <- LevN
                    sig[i] <- sig[i]/d
              }
       cutp[[i]] <- val$cuts
     }
  return(list(cutp=cutp,Disc.data=Disc))
}
