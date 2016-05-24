calc_N_Pot_Cuts_Sim <-
function(Daten, AnzCO){
  ld  <- ncol(Daten)
  cuts<-list()
  for(i in 2:ld){
      d<- sort(unique(round(Daten[,i], 1)))
      l<- length(d)
      
      if (l==1) cuts[[i]]<- d+0.1
      else{cuts[[i]]<-(d[1:l-1]+d[2:l])/2 
           if (l > AnzCO) cuts[[i]]<- sample(cuts[[i]], AnzCO)
    }}

  # Ausgabe:  
    return(cuts)
}
