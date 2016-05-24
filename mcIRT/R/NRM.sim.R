NRM.sim <-
function(Parlist,pp)

{
  
  # check input
  if(!is.list(Parlist)){stop("Parlist has to be a list!\n")}
  
  lparl <- sapply(Parlist,function(count) length(count) %% 2)
  if(any(lparl != 0)){stop(paste("Wrong number of categories for items:",which(lparl != 0)))}
  
  # start simulation
  
  realis <- lapply(Parlist,function(ITEM)
  {
    Pdist <- mapply(function(ze,la) exp(ze + la*pp),ze=ITEM[1:(length(ITEM)/2)],la=ITEM[(length(ITEM)/2 + 1):length(ITEM)])
    Pd2   <- Pdist / rowSums(Pdist)
    
    spalte <- sapply(1:length(pp),function(NEX)
    {
      #back <- sample(1:(length(ITEM)/2),1,prob=Pd2[NEX,])
      back <- sample(0:((length(ITEM)/2)-1),1,prob=Pd2[NEX,])
      back
    })
    #factor(spalte)
  })
  
  dfsim <- data.frame(realis)
  #dfsim <- as.matrix(realis)
  #
  dfsim
}
