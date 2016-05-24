#' Returns the country name corresponding to the passed iso code (3 letter, 2
#' letter or numeric).
#' 
#' Searches \code{getMap()@@data} to find the iso code. By default it returns the
#' string in the ADMIN column. By modifying nameColumn you can also get it to
#' return values from any other columns in \code{getMap()@@data} - see the examples.
#' Thus it can also be used to convert between ISO codes.
#' 
#' You could optionally provide a dataframe containing alternate iso
#' conversions using lookup= . The passe dataframe would need to contain at
#' least one of the following columns containing 2 letter, 3 letter or numeric
#' iso codes respectively : ISO_A2, ISO_A3, ISO_N3.
#' 
#' @param iso iso code to convert to a country name
#' @param lookup the dataframe containing iso codes and country names
#' @param nameColumn which column to get the name from, see examples
#' @return The country name (or other field) associated with the ISO code
#' passed.  NA is returned if no matching code is found.
#' @author Andy South
#' @keywords manip
#' @examples
#' 
#' isoToName('gb')
#' isoToName('gbr')
#' isoToName(826)
#' isoToName('uk') #generates a warning and returns NA
#' #beware that using nameColumn may be vulnerable to future changes 
#' #in column names in Natural Earth data
#' isoToName('gb',nameColumn='ABBREV') #returns abbreviation
#' isoToName('gb',nameColumn='ISO_A3') #returns iso3 for this iso2
#' isoToName('gbr',nameColumn='continent') #returns continent for this iso3
#' 
#' @export isoToName
`isoToName` <-
  
  function( iso = ""
          , lookup = getMap()@data  
          , nameColumn ='ADMIN'
            
          ){
    
    #iso can be a string of 2 or 3 chars or a number
    
    #get iso 2 or 3 from the string
    isoColumn <- NA
    if ( is.numeric(iso) )
        isoColumn <- 'ISO_N3'
    else if ( is.character(iso) && nchar(iso) == 3 )
        isoColumn <- 'ISO_A3'
    else if ( is.character(iso) && nchar(iso) == 2 )
        isoColumn <- 'ISO_A2' 
    else
    {
        warning(paste("iso= should be a 2 or 3 char string or a number, yours is",iso,"returning NA"))
        return(NA)
    }

    #if a character make sure it's in uppercase for comparison
    if ( is.character(iso) )  iso <- toupper(iso)
 
    index <- which( lookup[isoColumn] == iso )

    name <- NA
    if (length(index) == 0)
    {
       warning(paste("no matching name for",iso,"returning NA"))
    } else
    {  
       name <- as.character( lookup[,nameColumn][index] )
    }
    
    return(name)
  }

#testing
#isoToName('gb')
#isoToName('gbr')
#isoToName(826)
#isoToName('uk') #generates a warning and returns NA
#beware that using nameColumn may be vulnerable to future changes in column names in Natural Earth data
#isoToName('gb',nameColumn='ABBREV') #returns abbreviation
#isoToName('gb',nameColumn='ISO_A3') #returns iso3 for this iso2
#isoToName('gb',nameColumn='continent') #returns continent for this iso2
