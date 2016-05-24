 
sampleBy <- function(formula, frac=0.1, replace=FALSE, data=parent.frame(),
                     systematic=FALSE
                     ){
  ddd<-splitBy(formula, data=data)

  ddd<- lapply(ddd, function(dat){
    if (systematic==TRUE){
      idx <- seq(1,nrow(dat),1/frac)
    } else {  
      idx <- sort(sample(1:nrow(dat), size=round(frac*nrow(dat)), replace=replace))
    }
    dat[idx,]
  }) 
  ddd <-do.call("rbind",ddd)
  return(ddd)
}
