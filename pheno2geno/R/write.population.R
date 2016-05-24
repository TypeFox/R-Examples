#
# write.population.r
#
# Copyright (c) 2010-2013 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified Apr, 2013
# first written Apr, 2013
# Contains: write.population, writeSingleFile
#

#  write.population
#
# DESCRIPTION:
#  Writes geno/phenotypic files from R environment special object to disk.
# PARAMETERS:
#   - population - An object of class population .
#   - offspring - Core used to specify names of children phenotypic ("offspring_phenotypes.txt") and genotypic ("offspring_genotypes.txt") files.
#   - founders - Core used to specify names of founders phenotypic ("founders_phenotypes.txt") file.
#   - map - Core used to specify names of genetic ("map_genetic.txt") and physical ("map_physical.txt") map files.
#   - verbose - Be verbose
#   - debugMode - 1: Print our checks, 2: print additional time information
# OUTPUT:
#   Null
#
write.population <- function(population, offspring = "offspring", founders = "founders", map = "map", verbose = FALSE, debugMode = 0){

  ### checks
  check.population(population)

  ### file names
  fileFoundersPheno  <- paste(founders,"_phenotypes.txt",sep="")
  fileOffspringPheno <- paste(offspring,"_phenotypes.txt",sep="")
  fileOffspringGeno  <- paste(offspring,"_genotypes.txt",sep="")
  fileAnnotations    <- paste(offspring,"_annotations.txt",sep="")
  fileMapPhys        <- paste(map,"_physical.txt",sep="")
  fileMapGen         <- paste(map,"_genetic.txt",sep="")

  ### initializing  
  s <- proc.time()
  
  ### offspring phenotypes coukld be euther a matrix with phenotypes (that we can save) or a name of the file (then we just warn that the data is already in that file)
  ### if it is null - error!
  if(is.character(population$offspring$phenotypes)){ 
    cat("Phenotype data for offspring are already stored in:",population$offspring$phenotypes,". The function will not overwrite that file.")
  }else{ 
    writeSingleFile(population$offspring$phenotypes, "offspring phenotypes", fileOffspringPheno, errIfNotFound=TRUE, verbose=verbose) # REQUIRED
  }
  writeSingleFile(population$founders$phenotypes, "founder phenotypes", fileFoundersPheno, verbose=verbose)     # OPTIONAL
  writeSingleFile(population$annots, "annotation", fileAnnotations, verbose=verbose)                            # OPTIONAL
  writeSingleFile(population$offspring$genotypes, "offspring genotypes", fileOffspringGeno, verbose=verbose)    # OPTIONAL
  writeSingleFile(population$maps$physical, "physical map", fileMapPhys, verbose=verbose, col.names=FALSE)      # OPTIONAL
  writeSingleFile(population$maps$genetic, "genetic map", fileMapGen, verbose=verbose, col.names=FALSE)         # OPTIONAL

  #**********FINALIZING FUNCTION*************
  e <- proc.time()
  if(verbose && debugMode==2) cat("read.population finished after", (e-s)[3], "seconds.\n")
  invisible(NULL)
}

#  writeSingleFile
#
# DESCRIPTION:
#  Writes single geno/phenotypic file.
# PARAMETERS:
#   - dataMatrix - an object of class population
#   - filename - name of the file that will be processed
#   - verbose - be verbose
#   - ... - passed to write.table
# OUTPUT:
#   Null
#
writeSingleFile   <- function(dataMatrix, dataType, filename, errIfNotFound = FALSE, verbose=FALSE, ...){
  if(file.exists(filename)){
    cat("file: ",filename," already exists and will not be overwritten!\n")
   }else{
    if(missing(dataMatrix) || is.null(dim(dataMatrix))){
      if(errIfNotFound) stop("no",dataType,"found")
      return(cat("no",dataType,"found\n"))
    }
    write.table(dataMatrix,file=filename,sep="\t",quote=FALSE,...)
    if(verbose) cat(dataType,"saved in file:",filename,"\n")
  }
  invisible(NULL)
}
