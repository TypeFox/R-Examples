#gets ISO3 code from synonyms of a country name



#' Internal function for getting the ISO3 country code for a country name
#' synonymn.
#' 
#' Searches countrySynonyms to get the ISO3 code. If the name is not found NA
#' is returned. Allows joining of imperfect names to other country data in
#' joinCountryData2Map( joinCode='NAME' )
#' 
#' 
#' @param oddName country name that user wishes to find code for
#' @return the ISO3 code (3 letters) corresponding to the country name passed,
#' or NA if one is not found
#' @author Andy South
#' @references This was derived and used with permission from the Perl Locale
#' package.  \cr Locale::Codes::Country_Codes.  \cr Thanks to Sullivan Beck for
#' pulling this together.  \cr Data sources are acknowledged here : \cr
#' http://search.cpan.org/~sbeck/Locale-Codes-3.23/lib/Locale/Codes/Country.pod
#' @keywords manip
#' @examples
#' 
#' rwmGetISO3("vietnam")
#' 
#' @export rwmGetISO3
`rwmGetISO3` <- function( oddName ){
  
  oddName <- as.character( oddName )
  
  #first get the synonyms data
  data("countrySynonyms",envir=environment(),package="rworldmap")
  countrySynonyms <- get("countrySynonyms")
  
  #oldway using an input file
  #inFile <- "countrySynonyms.txt"
  #ncol <- max(count.fields(inFile, sep = "\t"))
  #countrySynonyms <- read.table(inFile,sep='\t', as.is=TRUE, fill=TRUE, header = FALSE
  #                              ,quote=""
  #                              ,col.names=c('ID','ISO3', paste("name", seq_len(ncol-2), sep = "")) )
   

  #firstNameColumn <- 2 #actually it's not 2 but won't do any harm to search through the ISO codes first
  
  #this works. indexing goes beyond end of rows for later columns
  #rowNum <- which( countrySynonyms[,2:length(countrySynonyms)] ==  oddName )
  
  #using tolower to avoid case things, start from col2 to include all synonyms
  #but tolower messes up dF into a char vector
  #rowNum <- which( tolower(dFcountries[,2:length(dFcountries)]) ==  tolower(oddName) )
  
  #keeping it simple
  cLow <- data.frame(lapply(countrySynonyms,tolower))
  
  rowNum <- which(cLow == tolower(oddName))
  
  #do I want it to return the oddName or NA if no match found ?
  
  #if no match no correct name
  if (length(rowNum)==0)
    {ISO3 <- NA} else 
    #this gets row nums for later columns      
    {
     #rowNum <- rowNum%%nrow(dFcountries)
     colNum <- 1+(rowNum-1)%/%nrow(countrySynonyms)   
     rowNum <- rowNum - (colNum-1)*nrow(countrySynonyms)    
     ISO3 <- countrySynonyms$ISO3[rowNum]
    }

  #I want as uppercase
  ISO3 <- toupper(ISO3)
  
  #just return the first one if there is more than one match
  return(ISO3[1])
}

#testing
#rwmGetISO3("Vietnam")
#rwmGetISO3("US")
#rwmGetISO3("USA")
#rwmGetISO3("Laos") 
#rwmGetISO3("Brunei") 
#rwmGetISO3("Republic of Lithuania")
#rwmGetISO3("Commonwealth of Dominica") #this is last one & can be a problem
#rwmGetISO3("China")
#rwmGetISO3("not a country") #returns what you passed to it. NA
