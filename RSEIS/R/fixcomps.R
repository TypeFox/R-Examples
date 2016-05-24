`fixcomps` <-
function(oldcomps, SEGY=FALSE)
  {
###################  compnent names
  ###################         in order to be consistent with the pickfiles
  ################      for matching, here we convert all compnent names to a standard convention
  ###########   the convention is: V N E    NOT Z
  ############  SEGY=456 usually  =>  4,5,6 = VNE
    #######       SEGY=123   =>  1,2,3  = VNE
  ####   sometimes acoustic A = 1, or infrasound I = 1
    
    if(missing(SEGY)) SEGY=FALSE

    UPcomps = toupper(oldcomps)
    comps = UPcomps
    nc = nchar(comps)
    
    comps =substr(comps, nc, nc)
    #####  convert any Z component to V temporarily for matching purposes
  
    comps[comps=="Z"] = "V"
    comps[comps=="U"] = "V"

    
    comps[UPcomps=="LD"] = "I"
   ###   comps[UPcomps=="L"] = "I"  ###  remove this option
  ###  comps[UPcomps=="A"] = "I"
    comps[UPcomps=="MIC"] = "I"

    
    
    ############  comps I  J K are for infrasound channels
    ########   I is for LARSON DAVIS?  by convention?
    

    
  if(SEGY==456)
    {
      #############  for SEGY data collected on refteks
      comps[comps=="4"] = "V" 
      comps[comps=="5"] = "N" 
      comps[comps=="6"] = "E" 
    }
  if(SEGY==123)
    {
      #############  for SEGY data collected on refteks
      comps[comps=="1"] = "V"
      comps[comps=="2"] = "N"
      comps[comps=="3"] = "E"
    }


    
  return(comps)

  }

