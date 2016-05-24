loadvfxmlbatch <- function( filename, patternMap, typeData = "vf" ) {
# loads CSV file with the filenames to find 
  csvforxml <- read.csv( filename, stringsAsFactors = FALSE )
  
# the CSV file must consist of two columns:
#  col 1: a suitable path for loading XML file, and
#  col 2: the type of subject "pwg" or "ctr"

# init
  xmlobject <- NULL
  xmlobject <- loadvfxml( csvforxml[1,1], patternMap = patternMap, typeData = typeData, typeSubject = csvforxml[1,2], extractionType = c( "average" ) )
  if( nrow( csvforxml ) == 1 ) return( xmlobject )
  for( i in 2:nrow( csvforxml ) ) {
    xmlobject[i,] <- loadvfxml( csvforxml[i,1], patternMap = patternMap, typeData = typeData, typeSubject = csvforxml[i,2], extractionType = c( "average" ) )
  }

  return( xmlobject )

}
  