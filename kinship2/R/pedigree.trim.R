# Automatically generated from all.nw using noweb

pedigree.trim <- function(removeID, ped){
  ## trim subjects from a pedigree who match the removeID 
  ## trim relation matrix as well

  if(class(ped) != "pedigree")
    stop("Must be a pegigree object.\n")

  rmidx <- match(removeID, ped$id)
  if(length(rmidx)>0) {
    pedtrimmed <- ped[-rmidx]
    return(pedtrimmed)
  } else {
    return(ped)
  }
}
 
