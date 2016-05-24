# Copyright
# This function is from the PGE Toolbox (Matlab)
# slightly modified
#Cai JJ (2008) PGEToolbox: A Matlab toolbox for population genetics and evolution
#Journal of Heredity Jul-Aug;99(4):438-40. doi:10.1093/jhered/esm127
#modified

hudsonkaplan85rm <- function(bial,populations) {

npops    <- length(populations)
sg.sites <- get_segsites(bial,populations)
RM       <- numeric(npops)

for(xx in 1:npops){
  
  m      <- length(sg.sites[[xx]])
  if(m<2){next}
  s.bial <- bial[,sg.sites[[xx]],drop=FALSE]
  Dcount <- matrix(0,m,m)
  pairs  <- combn(m,2)
  
  for(yy in 1:dim(pairs)[2]){

     id1     <- pairs[,yy][1]
     id2     <- pairs[,yy][2]
     site1_2 <- s.bial[,c(id1,id2),drop=FALSE]
     numHap  <- dim(unique(site1_2))[1]
     Dcount[id1,id2] <- numHap
  }
     koord   <- which(Dcount==4,arr.ind=TRUE)
     if(length(koord)==0){next}
     x       <- koord[,1]
     y       <- koord[,2]
   
  
    idx <- vector(,length(x)+1)
    x<-c(x,0); y<-c(y,0); idx[length(idx)]<-TRUE
    while (any(idx)){
        id  <- which(idx)
        x   <- x[-id]
        y   <- y[-id] 
        idx <- idx[-id]
    if (length(x)>=2){
     for (k in 2:length(x)){
         if(x[k]>=x[k-1] && y[k]<=y[k-1]){
            idx[k-1]<-TRUE
         }
     }
    }
    }

    if (length(x)>1){
        #x=[x(end:-1:1);0]; y=[y(end:-1:1);0]; idx(end)=1;
      x <- c(x[length(x):1],0); y <- c(y[length(y):1],0); idx <- c(idx,TRUE)
    while (any(idx)){
        id  <- which(idx)
        x   <- x[-id]
        y   <- y[-id] 
        idx <- idx[-id]
    if (length(x)>=2){
     for (k in 2:length(x)){
         if(x[k]>=x[k-1] && y[k]<=y[k-1]){
            idx[k-1]<-TRUE
         }
     }
    }
    }          
    }

RM[xx] <- length(x)

}

return(RM)
}
