NLM.sim <-
function(Parlist,pp)

{
  # CONTROL Parlist
  realis <- lapply(Parlist,function(ITEM)
  {
    # prob of solution
    Ploes     <- 1/(1+exp(-ITEM[2] - ITEM[1]*pp)) 
    distrpar  <- ITEM[-(1:2)]
    distrparL <- split(distrpar,rep(1:2,each=length(distrpar)/2))
    Pdist <- mapply(function(ze,la) exp(ze + la*pp),ze=distrparL[[1]],la=distrparL[[2]])
    Pd2   <- Pdist / rowSums(Pdist)
    
    spalte <- sapply(1:length(pp),function(NEX)
    {
      rf_real <- sample(c(1,0),1,prob=c(Ploes[NEX],1-Ploes[NEX]))
      
      if(rf_real == 0)
      {
        back <- sample(2:(ncol(Pd2)+1),1,prob=Pd2[NEX,]) 
        back - 1 # um auf dieselbe form zu kommen wie bei nrm
      } else 
      {
        rf_real  - 1
      }
    })
    #factor(spalte)
  })
  
  data.frame(realis)
}
