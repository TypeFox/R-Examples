ComputeMarketPrices<-function(Availability,OffresPrix,Conso,OffresPrixb=NULL)
{
  stopifnot(sum(is.na(Conso))==0)
  stopifnot(sum(Availability<0)==0)
  stopifnot(dim(Availability)==dim(OffresPrix))
  stopifnot(dim(Availability)[1]==length(Conso))
  if (is.null(OffresPrixb))
  {
    ####
    #### solves for all time step t
    #### min P_{t,x} sum_x P_{t,x}OffresPrix[t,x]
    #### sum_x P_{t,x}=Conso[t]
    #### 0<=P_{t,x}<=availability[t,x]
    
    if (is.vector(Conso))
    {
      res=OptimPriceMarket_l(OffresPrix,Availability,Conso)
    }else{
      ## not implemented xts implementation ? 
      stop('Conso should be a vector')
    }
    
  }else
  {
    ####
    #### solves for all time step t
    #### min P_{t,x} sum_x P_{t,x}(P_{t,x}OffresPrix[t,x]+OffresPrixb[t,x])
    #### sum_x P_{t,x}=Conso[t]
    #### 0<=P_{t,x}<=availability[t,x]
    
    if (is.vector(Conso))
    {
      res=OptimPriceMarket_q(OffresPrix,OffresPrixb,Availability,Conso)
    }else{
      ## not implemented xts implementation ? 
      stop('Conso should be a vector')
    }
  }
return(res)
}





OptimPriceStorage<-function(Prices,Pplus,Pmoins,Cplus,Cmoins=0,
                            efficiencyS=0,efficiencyP=efficiencyS,networkTax=0)
{
  stopifnot(sum(is.na(Prices))==0)
  stopifnot(sum(Prices<0)==0)
  if (is.vector(Prices))
  {
    if (length(Pplus)==1) Pplus=array(Pplus,length(Prices))
    if (length(Pmoins)==1) Pmoins=array(Pmoins,length(Prices))
    if (length(Cplus)==1) Cplus=array(Cplus,length(Prices))
    if (length(Cmoins)==1) Cmoins=array(Cmoins,length(Prices))
    if (efficiencyS==0&efficiencyP==0&networkTax==0)
    {
      Operation=OptimPriceStorage_(Prices,Pmoins,Pplus,Cmoins,Cplus)$xEtoile
      return(list(Operation=Operation,
                  Revenue=-Prices*Operation))
    }else
    {
      CCPWLfuncList=new(cplfunctionvec)
      CCPWLfuncList$SerialPush_2Breaks_Functions(Prices*efficiencyP,
                                                 Prices/efficiencyS+networkTax,
                                                 array(-Inf,length(Cmoins)),
                                                 array(0,length(Cmoins)));
      Operation=CCPWLfuncList$OptimMargInt(Pmoins,Pplus,Cmoins,Cplus)$xEtoile
      return(list(Operation=Operation,
                  Revenue=-(Prices[Operation<0]*efficiencyP*Operation[Operation<0]+(Prices[Operation>0]/efficiencyS+networkTax)*Operation[Operation>0])))
    }
    
    
  }else{
    ## not implemented xts implementation ? 
    stop('Prices should be a vector')
  }
  
}
