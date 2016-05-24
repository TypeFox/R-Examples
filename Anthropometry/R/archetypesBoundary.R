archetypesBoundary <- function(data,numArch,verbose,numRep){

  ldata <- data 
 
  #Run archetypes algorithm repeatedly from 1 to numArchet archetypes:
  sequen <- seq(length = numArch)
  lass <- stepArchetypesRawData(data = ldata, numArch = sequen, 
                            numRep = numRep, verbose = verbose) 

  return(lass) 
}











