codes2names <- function(x, variable="circumstances", type="euring"){
    # x: variable with codes as given in the EURING data base
    # variable: character, one of "circumstances", "conditions", "schemes", "species"
    # type: charcter, one of "euring" or "bto" (only in use for circumstances)
  
  #filename <- system.file("data", paste0("codes_", variable,".txt"), package = "birdring")
  #data(circumstances, envir = environment())
  #data(conditions, envir = environment())
  #data(schemes, envir = environment())
  #data(species, envir = environment())
  type <- tolower(type)  # make sure that all letters are lower letters
  codesdat <- get(variable) #read.table(filename, header=TRUE)
  if(variable=="schemes") {
    codesdat$Name <- paste(codesdat$Country, codesdat$Centre)
    codesdat$Code <- gsub(" ", "", codesdat$Code)
  }
  if(variable=="conditions") codesdat$Name <- codesdat$Description
  
  if(variable=="circumstances"&type=="bto"){
    codesdat$Name <- as.character(codesdat$Name)
    codesdat$Name[codesdat$Code==20] <- "caught by ringer" 
    codesdat$Name[codesdat$Code==27] <- "caught in nestbox" 
    codesdat$Name[codesdat$Code>=10&codesdat$Code<=19] <- "shot" 
    codesdat$Name[codesdat$Code>=20&codesdat$Code<=26] <- "trapped" 
    codesdat$Name[codesdat$Code==28] <- "ring read in field" 
    codesdat$Name[codesdat$Code==29] <- "colour marks seen" 
    codesdat$Name[codesdat$Code==30] <- "oiled" 
    codesdat$Name[codesdat$Code==31] <- "pollution" 
    codesdat$Name[codesdat$Code==32] <- "on wire or netting" 
    codesdat$Name[is.element(codesdat$Code, c(33, 34))] <- "in net or cage" 
    codesdat$Name[codesdat$Code==35] <- "electrocuted" 
    codesdat$Name[is.element(codesdat$Code, c(36, 37,38))] <- "poisoned" 
    codesdat$Name[codesdat$Code==40] <- "hit by car" 
    codesdat$Name[codesdat$Code==41] <- "hit by train" 
    codesdat$Name[codesdat$Code==42] <- "hit by plain" 
    codesdat$Name[codesdat$Code==43] <- "hit wires" 
    codesdat$Name[codesdat$Code==44] <- "hit glass" 
    codesdat$Name[codesdat$Code==45] <- "hit building" 
    codesdat$Name[codesdat$Code==46] <- "in building" 
    codesdat$Name[is.element(codesdat$Code, c(47, 48))] <- "accidental" 
    codesdat$Name[codesdat$Code==49] <- "drowned" 
    codesdat$Name[is.element(codesdat$Code, c(50, 51))] <- "injury" 
    codesdat$Name[is.element(codesdat$Code, c(52:59))] <- "disease" 
    codesdat$Name[is.element(codesdat$Code, c(60, 66:69))] <- "predated" 
    codesdat$Name[codesdat$Code==61] <- "cat" 
    codesdat$Name[codesdat$Code==62] <- "domestic animal" 
    codesdat$Name[codesdat$Code==63] <- "wild animal" 
    codesdat$Name[is.element(codesdat$Code, c(64,65))] <- "bird of prey" 
    codesdat$Name[codesdat$Code==70] <- "drowned" 
    codesdat$Name[is.element(codesdat$Code, c(71, 72, 73, 75, 76, 79))] <- 
      "natural causes" 
    codesdat$Name[codesdat$Code==74] <- "cold weather" 
    codesdat$Name[codesdat$Code==77] <- "ice" 
    codesdat$Name[codesdat$Code==78] <- "storm" 
    codesdat$Name[codesdat$Code==81] <- "colour rings seen" 
    codesdat$Name[codesdat$Code==82] <- "neck collar seen" 
    codesdat$Name[codesdat$Code==83] <- "wing-tag seen" 
    codesdat$Name[codesdat$Code==84] <- "radio-tagged" 
    codesdat$Name[codesdat$Code==85] <- "satellite-tagged" 
    codesdat$Name[codesdat$Code==86] <- "transponder tag" 
    codesdat$Name[codesdat$Code==2] <- "ring only" 
    codesdat$Name[codesdat$Code==3] <- "leg only" 
    codesdat$Name[codesdat$Code==6] <- "on ship" 
  }
  
  namesx <- codesdat$Name[match(as.character(x), as.character(codesdat$Code))][drop=TRUE]
  return(namesx)
}