extendChi2 <-
function(data,alp=0.5){
        #Phase 1
        dat <- data
        p <- ncol(dat)-1      
        alpha <- alp; d=0.1

        XiO <- Xi(dat)
        Disc <- chiM(dat,alpha=alpha)$Disc.data
        XiN <- Xi(Disc)
        e <- 0.01
     while((XiN-XiO)<e){
               if(alpha<=0.2) break
               alpha <- alpha-d
               Disc <- chiM(Disc,alpha=alpha)$Disc.data
               XiN <- Xi(Disc) 
               XiO <- XiN               
    }

#Phase 2
       options(digits=3)
       eps <- 0.0001
       sig <- array(alpha,p)
       d <- 10
    cutp=list()
    for(i in 1:p){
                XiO <- Xi(Disc)
              while(TRUE){
                    val <- value(i,Disc,alpha=sig[i])
                    Disc <- val$disc
                    XiN <- Xi(Disc)#inconRate for data
                    XiO <- XiN
                    sig[i] <- sig[i]/d
                    if(XiN<XiO|| sig[i]<=eps) break 
              }
      cutp[[i]]<- val$cuts
     }
  return(list(cutp=cutp,Disc.data=Disc))
}
