#' Produce regional data from country level data
#' 
#' A function to aggregate country level data into regional data. For example
#' finding the total population of Asia, Europe,etc, from country level
#' populations. As well as sums, other functions can be used, like mean,
#' median, min, max, etc. There are currently 8 choices of region and 4 choices
#' of country code.
#' 
#' The user must specify 'nameJoinColumn' from their data which contains
#' country codes, and joinCode which specifies the type of code. regionType
#' specifies which regions to aggregate the data to. Using FUN='identity' will
#' return the neames of the countries within each region.
#' 
#' @param regionType Must be one of: "GEO3", "GEO3major", "IMAGE24", "GLOCAF",
#' "Stern", "SRES", "SRESmajor" or "GBD"
#' @param inFile a data frame
#' @param nameDataColumn The name of the data column to aggregate
#' @param joinCode The type of code to join with. Must be one of: "ISO2",
#' "ISO3", "Numeric" or "FIPS"
#' @param nameJoinColumn The name of a column of inFile. Contains joining
#' codes.
#' @param FUN A function to apply to each region, e.g. 'mean'
#' @param \dots further arguments to be passed to FUN, e.g. na.rm=TRUE
#' @return If FUN returns a single value, country2Region returns a data frame,
#' with value of FUN for each region.
#' 
#' If FUN returns more than one value, country2Region will return a list, with
#' one element for each region.
#' @seealso For producing maps of regional data from aggregated country level
#' data, see \code{\link{mapByRegion}}
#' @keywords manip
#' @examples
#' 
#' data(countryExData)
#' 
#' #to report which countries make up regions
#' country2Region(regionType="Stern")
#' 
#' #Using country2Region to calculate mean Environmental Health index in Stern regions.
#' sternEnvHealth <- country2Region(inFile=countryExData
#' 		,nameDataColumn="ENVHEALTH"
#' 		,joinCode="ISO3"
#' 		,nameJoinColumn="ISO3V10"
#' 		,regionType="Stern"
#' 		,FUN='mean'
#' 		)
#' 
#' print(sternEnvHealth)
#' 
#' #A simple plot of this data.
#' #dotchart(sort(sternEnvHealth))
#' dotchart(sort(sternEnvHealth[,1]))
#' 
#' #use FUN='identity' to see which countries in your data belong to which region.
#' country2Region(inFile=countryExData
#' 		,nameDataColumn="Country"
#' 		,joinCode="ISO3"
#' 		,nameJoinColumn="ISO3V10"
#' 		,regionType="Stern"
#' 		,FUN='identity'
#' 		)
#' 
#' #Change FUN to length, to count the number of countries in each region.
#' country2Region(inFile=countryExData
#' 		,nameDataColumn="Country"
#' 		,joinCode="ISO3"
#' 		,nameJoinColumn="ISO3V10"
#' 		,regionType="Stern"
#' 		,FUN='length'
#' 		)
#' 
#' 
#' 
#' @export country2Region
`country2Region` <-
function(regionType = ''
        ,inFile = ''
        ,nameDataColumn = ''
        ,joinCode = ''
        ,nameJoinColumn = ''
        ,FUN=mean
        ,...)
{
#the data specifiying which countries are in which regions
data("countryRegions",envir=environment(),package="rworldmap")
countryRegions <- get("countryRegions")

valid_classification_types<-c("REGION","continent","GEO3","GEO3major","IMAGE24","GLOCAF","Stern","SRESmajor","SRES","GBD","AVOIDnumeric","AVOIDname","LDC","SID","LLDC")

#prompt the user for a regionType if one is not specified    
while(!(regionType %in% valid_classification_types))
   {
    regionTypeOptions <- paste(valid_classification_types,collapse=" ")
    regionType <- readline(paste("Please enter a valid regionType. The options are:\n",regionTypeOptions,"\n"))
   }
    

#if only a regionType is specified then return which countries belong to it
if ( !is.data.frame(inFile) && inFile == '' && nameDataColumn == '' && joinCode == '' && nameJoinColumn == '' )
   {
    message(paste('Countries are classified into the following ', regionType, 'regions'))
    FUN <- 'identity'
    tapply(countryRegions$ADMIN,countryRegions[[regionType]],FUN,...)  
    
   } else
   {
    #copying the dataFrame passed to the function
    dF<-inFile  
    
    #checking that the function is a character string
    if( !is(FUN,"character") )
      {
       warning(paste("option FUN should be a character string enclosed in quotes, e.g. 'mean', 'sum' using 'mean'") )
       FUN = "mean"
      }
    
    #valid_code_types<-c("ISO3","ISO2","Numeric","Name","FIPS")
    valid_code_types<-c("ISO3","ADMIN")
    
    #checking that the country join code is valid
    if(!(joinCode %in% valid_code_types))stop("joinCode is invalid. The options are: ",paste(valid_code_types,collapse=" "))
    
    #check whether nameJoinColumn is in the user data
    ## check that the column name exists in the data frame
    if ( is.na(match(nameDataColumn, names(dF)) )){
      stop("your chosen nameDataColumn :'",nameDataColumn,"' seems not to exist in your data, columns = ", paste(names(dF),""))
      return(FALSE)
    }     
    if ( is.na(match(nameJoinColumn, names(dF)) )){
      stop("your chosen nameJoinColumn :'",nameJoinColumn,"' seems not to exist in your data, columns = ", paste(names(dF),""))
      return(FALSE)
    } 
   
    #Create a temporary, simple lookup table from the master look up table.
    #This just contains the country codes and groupings asked for in the 1st and 2nd column respectively.
    subLookUpTable<-countryRegions[,c(joinCode,regionType)]
    
    #Countries will not always have a code for every code type.
    #e.g. Palestine has a code under ISO3 but not FIPS, which has seperate codes for the West Bank and the Gaza Strip.
    #This removes rows that have an NA in the code column. Otherwise NAs in the data column will match with it,
    #causing mis-classification.
    subLookUpTable<-subLookUpTable[!is.na(subLookUpTable[,joinCode]),]
    
    classified<-subLookUpTable[match(dF[,nameJoinColumn],subLookUpTable[,joinCode]),regionType]
    
    output <- tapply(dF[,nameDataColumn],classified,FUN,...)
    
    #looking at returning in a more useful format
    idString <- paste(FUN,nameDataColumn,'by',regionType,sep='')
    dFout <- data.frame( x=output )
    names(dFout)[1] <- idString
    
    #if there is more than one element per region (e.g.identity) return as a list
    #else return as a dat frame, easier to access & output
    if ( is.list(dFout[,1]) )
         return(dFout[,1])
    else return(dFout)    
    
    #write.csv(output,'test.csv') #to output to a csv file
    
    
    
   } #end of if inFile etc. are specified

} #end of country2Region

