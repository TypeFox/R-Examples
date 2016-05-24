ConvertConcentration <- function(x, convert, percent) {

  ## define conversion functions
  
  DryToWet <- function(dry, percent.moisture) {

    wet <- dry * (1 - (percent.moisture / 100))

    return(wet)

  }

  WetToDry <- function(wet, percent.moisture) {

    dry <- wet / (1 - (percent.moisture / 100))

    return(dry)

  }

  WetToLipid <- function(wet, percent.lipid) {

    lipid <- wet * (100 / percent.lipid)

    return(lipid)

  }

  LipidToWet <- function(lipid, percent.lipid) {

    wet <- lipid * (percent.lipid / 100)

    return(wet)

  }

  ## error check

  if(!(convert == "dry.to.wet" | convert == "wet.to.dry" |
       convert == "wet.to.lipid" | convert == "lipid.to.wet")) {
    
    stop("the convert argument is not vaild")
    
  }

  ## make calculation

  if(convert == "dry.to.wet") { y <- DryToWet(x, percent) }

  if(convert == "wet.to.dry") { y <- WetToDry(x, percent) }

  if(convert == "wet.to.lipid") { y <- WetToLipid(x, percent) }

  if(convert == "lipid.to.wet") { y <- LipidToWet(x, percent) }

  return(y)

}

    

    

    
