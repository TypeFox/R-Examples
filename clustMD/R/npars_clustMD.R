npars_clustMD <-
function(model, D, G, J, OrdIndx){
    if(model=="EII"){
      npars <- 1 + (G-1) + G*D
      # if nominal
      if(J > OrdIndx){
        npars <- 1 + (G-1) + G*OrdIndx + (G-1)*(D-OrdIndx)
      }
    }else if(model=="VII"){
      npars <- G + (G-1) + G*D
      # if nominal
      if(J > OrdIndx){
        npars <- G + (G-1) + (G-1) + G*OrdIndx + (G-1)*(D-OrdIndx)
      }
    }else if(model=="EEI"){
      npars <- 1 + (D-1) + (G-1) + G*D 
      # if nominal
      if(J > OrdIndx){
        npars <- 1 + (OrdIndx-1) + (G-1) + G*OrdIndx + (G-1)*(D-OrdIndx)
      }
    }else if(model=="VEI"){
      npars <- G + (D-1) + (G-1) + G*D 
      # if nominal
      if(J > OrdIndx){
        npars <- G + (G-1) + (OrdIndx - 1) + (G-1) + G*OrdIndx + (G-1)*(D-OrdIndx)
      }
    }else if(model=="EVI"){
      npars <- 1 + G*(D-1) + (G-1) + G*D
      # if nominal
      if(J > OrdIndx){
        npars <- 1  + G*(OrdIndx - 1) + (G-1)*(D - OrdIndx - 1) + (G-1) + G*OrdIndx + (G-1)*(D-OrdIndx)
      }
    }else if(model=="VVI"){
      npars <- G + G*(D-1) + (G-1) + G*D
      # if nominal
      if(J > OrdIndx){
        npars <- G + (G-1) + G*(OrdIndx-1) + (G-1)*(D - OrdIndx - 1) + (G-1) + G*OrdIndx + (G-1)*(D-OrdIndx)
      }
    }
    npars
  }
