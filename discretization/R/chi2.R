chi2 <-
function(data,alp=0.5,del=0.05){
#Phase 1
    dat <- data
    p <- ncol(dat)-1         
    
    alpha <- alp; delta <- del; d <- 0.1 
    inconRate <- incon(dat)
    Disc <- chiM(dat,alpha=alpha)$Disc.data
          
       while(inconRate<delta){
           if(alpha<=0.2) break
           alpha <- alpha-d
           Disc <- chiM(Disc,alpha=alpha)$Disc.data
               inconRate <- incon(Disc) 
       }
       
#Phase 2
    options(digits=3)
    eps <- 0.01
    sig <- array(alpha,p)
    d <- 0.05
    cutp=list()
    for(i in 1:p){
        while(TRUE){
             val <- value(i,Disc,alpha=sig[i])
             Disc <- val$disc
             inconRate <- incon(Disc)#inconRate for data
             sig[i] <- sig[i]-d
             if(inconRate>delta || sig[i]<=eps) break 
        }
     cutp[[i]] <- val$cuts
    }
return(list(cutp=cutp,Disc.data=Disc))
}
