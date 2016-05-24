expandCensus<-function(popdata,year.range=range(as.integer(names(popdata))),years=as.integer(names(popdata))){


  times <- c(year.range[1], sort(years), year.range[2])
  times <- as.numeric(times)
  inter <- diff(times)/2
  nseq <- 1:length(inter) - 1
  mseq <- 2:length(inter)
  interval <- inter[mseq] + inter[nseq]
  interval[1] <-interval[1] + inter[1]
  interval[length(interval)] <-interval[length(interval)] + inter[length(inter)]
  interval<-ceiling(interval)

  newp<-list()
  index=1
  for (k in interval){
  
    if(index==1) n<-1 else n <- length(newp)+1
    for (i in n:(n+k-1)){
      newp[[i]]<-popdata[[index]]
    }
  index=index+1
 }
 
  names(newp) = year.range[1]:year.range[2]


newp
}