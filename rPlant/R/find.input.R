find.input <- function(x){
  if (x == "bed"){
    return("inputBED")
  } else if (x == "bim"){
    return("inputBIM")
  } else if (x == "fam"){
    return("inputFAM")
  } else if (x == "tfam"){
    return("inputTFAM")
  } else if (x == "tped"){
    return("inputTPED")
  } else if (x == "map"){
    return("inputMAP")
  } else if (x == "ped"){
    return("inputPED")
  } else {
    return("Error: wrong input files")
  }
}
