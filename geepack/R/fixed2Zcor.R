##
## Generate zcor vector from 
## 1) fixed correlation matrix
## 2) id information
## 3) waves information

## The zcor-vector contrains entries only for clusters 
## of size larger than 1 

fixed2Zcor <- function(cor.fixed, id, waves){
  zcor<-NULL
  cnt <- 1
  uniq.id <- unique(id)
  for (ii in uniq.id){
    cwaves <- waves[id==ii]
    if (length(cwaves)>1) {
      for (kk in 1: (length(cwaves)-1)) {
        for (mm in (kk+1) : length(cwaves))  {
            vvv <- cor.fixed[cwaves[mm],cwaves[kk]]
           zcor<-c(zcor,vvv)

        }
      }
    }
  }
  zcor
}