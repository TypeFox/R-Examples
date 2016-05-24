#' Importing a Genepop file
#' 
#' @description This function converts a Genepop file into an object 
#' with a genetic matrix (G) and a structures matrix (S).
#' 
#' @param genefile Genepop file.
#' 
#' @return A list with the objects G (genetic matrix) and S (structures matrix).
#'
#' @examples 
#' \dontrun{
#' # ingpop, file with Genepop format in the folder "/extdata" of the package
#' 
#' ecopath <- paste(path.package("EcoGenetics"), "/extdata/ingpop", sep = "")
#' ingpop <- eco.genepop2df(ecopath)
#' ingpop
#' }
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}, adapting code
#' written by Emiel van Loon and Scott Davis
#' 
#' @export


eco.genepop2df <- function(genefile = NULL) {

  # read the file into memory
  con <- file(description = genefile, open = "rb")
  lines <- readLines(con)
  close(con)
  fileLength = length(lines)
  
  # checking that the file ends with a newline
  endline <- lines[fileLength]
  endline <- strsplit(endline, " ")[[1]]
  if(any(endline != "")) {
    stop("the file does not end with a newline")
  }
  
  #pop locations ignoring first 'title' line
  popLocations <- grep("POP", lines[2:fileLength], ignore.case=TRUE)
  popLocations <- popLocations + 1 
  
  # get title from first line
  # title <- lines[1]
  
  #get Loci Column names
  lociNames <- lines[2:(popLocations[1]-1)]
  lociNames <- gsub("\t","",lociNames)
  
  # number of individuals and populations
  npop <- length(popLocations)
  initial <- popLocations[1] 
  end <- length(lines)
  rownumber <- end - initial - npop
  
  nombres.loci <- unlist(strsplit( lociNames, '[[:space:]]+'))
  colnumber <- length(nombres.loci)
  
  out <- matrix(0, ncol = colnumber, nrow = rownumber)
  estructuras <- rep(0, rownumber)
  nombres <- character()
  
  #parsing of data between pops into genepop
  counter <- 1
  
  for(i in 1:npop) {
    beginLine <- popLocations[i] + 1
    endLine <- 0
    
    if( i == npop) {  
      endLine <- length(lines) - 1
      
    } else {	
      endLine <- popLocations[i+1] - 1
    }
 
    #parse individual line
    
    for(line in lines[beginLine:endLine]) {
      
      #split id & alleles apart
      individuo <- unlist(strsplit(line, ","))
      nombres[counter]<- individuo[1]
      estructuras[counter] <- i
      haplotipo <- individuo[2]
      
      #split alleles apart on whitespace
      loci <- unlist(strsplit(haplotipo, '[[:space:]]+'))
      loci <- loci[2:length(loci)]
      
      #saving non zero lines into out
      
      out[counter, ] <- loci
      counter <- counter + 1 
    }
  }
  out <- data.frame(out)
  rownames(out) <- nombres
  colnames(out) <- nombres.loci
  estructuras <- data.frame(factor(estructuras))
  colnames(estructuras)[1] <- "pop"
  rownames(estructuras) <- nombres
  list("G" = out, "S" = estructuras)
}
