#' Joins user country referenced data to a map
#' 
#' Joins user data referenced by country codes or names to an internal map,
#' ready for plotting using \code{\link{mapCountryData}}.  Reports join
#' successes and failures.
#' 
#' Joins data referenced by country codes to an internally stored map to enable
#' plotting.  The user specifies which country code their data are referenced
#' by, and the name of the column in their data containing that referencing
#' data. The user can choose from different map resolutions, using the function
#' \code{\link{getMap}} to retrieve the map. The function reports on how many
#' countries successfully join to the map. Data can then be plotted using
#' \code{\link{mapCountryData}}. NEW to version 1.01 Oct 2012 : for
#' joinCode='NAME' alternative country names are matched using
#' \code{\link{countrySynonyms}}.
#' 
#' The projection argument has now been deprecated, you can project maps using
#' package rgdal as shown below and in the FAQ. \cr library(rgdal) \cr #first
#' get countries excluding Antarctica which crashes spTransform \cr sPDF <-
#' getMap()[-which(getMap()$ADMIN=='Antarctica'),] \cr #transform to robin for
#' the Robinson projection \cr sPDF <- spTransform(sPDF, CRS=CRS("+proj=robin
#' +ellps=WGS84")) \cr mapCountryData( sPDF, nameColumnToPlot="REGION")
#' 
#' @param dF R data frame with at least one column for country reference and
#' one column of data
#' @param joinCode how countries are referenced options
#' "ISO2","ISO3","FIPS","NAME", "UN" = numeric codes
#' @param nameJoinColumn name of column containing country referencing
#' @param nameCountryColumn optional name of column containing country names
#' (used in reporting of success/failure)
#' @param suggestForFailedCodes NOT YET ENABLED T/F whether you want system to
#' suggest for failed codes
#' @param projection DEPRECATED JUNE 2012
#' @param mapResolution resolution of the borders in the internal map, only for
#' projection='none' : options 'low', 'medium'
#' @param verbose if set to FALSE it doesn't print progress messages to console
#' @return An R 'SpatialPolygonsDataFrame' [package "sp"] object with the
#' passed data joined to it
#' @author andy south
#' @seealso \code{\link{mapCountryData}}, \code{\link{getMap}}
#' @keywords dplot
#' @examples
#' 
#' data("countryExData",envir=environment(),package="rworldmap")
#' 
#' sPDF <- joinCountryData2Map(countryExData
#'               , joinCode = "ISO3"
#'               , nameJoinColumn = "ISO3V10"
#'               )
#' mapCountryData( sPDF
#'               , nameColumnToPlot="BIODIVERSITY" 
#'               )
#' 
#' 
#' @export joinCountryData2Map
`joinCountryData2Map` <-
function( dF
        , joinCode = "ISO3" #options "ISO2","ISO3","FIPS","NAME","UN"
        , nameJoinColumn = "ISO3V10"
        , nameCountryColumn = "Country"
        , suggestForFailedCodes = FALSE 
        , mapResolution="coarse" 
        , projection=NA 
        , verbose = FALSE #if set to FALSE it doesn't print progress messages to console                         
        )
   {

    mapWithData <- getMap(resolution=mapResolution)

		if ( ! is.na(projection))
		   warning('the projection argument has been deprecated, returning Lat Lon, use spTransform from package rgdal as shown in help details or the FAQ')
    
    #test whether user joinCode is one of permitted
    #natural earth data has : "ISO_A2","ISO_A3","FIPS_10_","ADMIN","ISO_N3",  #*beware that "NAME" in nat earth has abbreviations that don't match  
    #here I allow Nat Earth codes or my original codes
    
    listJoinCodesNew <- c("ISO_A2","ISO_A3","FIPS_10_","ADMIN","ISO_N3" )
    listJoinCodesOld <- c("ISO2","ISO3","FIPS","NAME","UN" )
    listJoinCodes <- c(listJoinCodesOld,listJoinCodesNew)
    if (joinCode %in% listJoinCodes == FALSE)
       {
        stop("your joinCode (",joinCode,") in joinCountryData2Map() is not one of those supported. Options are :",paste(listJoinCodes,""),"\n")
        return(FALSE)
       }
    
    #converting any join codes from old to new
    joinCodeOld <- joinCode
    if (joinCode %in% listJoinCodesOld)
       {joinCode <- listJoinCodesNew[match(joinCode,listJoinCodesOld)]}
      
    
    ## check that the join column exists in the user data
    if ( is.na(match(nameJoinColumn, names(dF)) )){
      stop("your chosen nameJoinColumn :'",nameJoinColumn,"' seems not to exist in your data, columns = ", paste(names(dF),""))
      return(FALSE)
    }       
       

    #dF2 <- merge.data.frame(dF, dFlookupCodes, by=nameJoinColumn)
    
    #using match rather than merge, faster and enables greater reporting of success & failure
    
    #match returns a vector of the positions of (first) matches of its first argument in its second. 
    #so perhaps I would also want to check that codes aren't repeated no can't do everything for people
    #!also want to find a way of coping with Namibia, the code NA gets interpreted as no data
           
    #copy the users nameJoinColumn to a new column named the same as the column in the map for the join code
    #e.g if user has ISO3166_3 it will be copied to ISO3
    dF[[joinCode]] <- as.character(dF[[nameJoinColumn]])
    
    #this removes any trailing spaces from the user data which could be a problem
    dF[[joinCode]] <- gsub("[[:space:]]*$","",dF[[joinCode]])
    
    #[:space:] is a pre-defined character class that matches space characters in your locale. 
    #* says to repeat the match zero or more times and $ says to match the end of the string.
    
    #23/5/12 if using NAME I could convert to ISO3 first using synonyms and match based on that
    #and set nameCountryColumn to what the join column was
    #but does everything become a bit hidden then ?
    #not really keeps it fairly simple
    if (joinCode=='ADMIN' ) #now it is ADMIN because I converted above from NAME to ADMIN
        { 
         #get the equivalent ISO 3 codes for the column
         #create new column
         dF$ISO3 <- NA
         for(i in 1:nrow(dF)) dF$ISO3[i] = rwmGetISO3( dF[[joinCode]][i] )
         #set join code to ISO3
         joinCode='ISO3';
         #set name for the country column to what user had as the name join column 
         nameCountryColumn=nameJoinColumn; 
         #and set the nameJoinColumn to ISO3 which we've added to the user data
         #nameJoinColumn='ISO3'
         }
    
    matchPosnsInLookup <- match(as.character(dF[[joinCode]])
                              , as.character(mapWithData@data[[joinCode]]))


    #count the NAs to find user countries that have failed to match
    failedCodes <- dF[[joinCode]][is.na(matchPosnsInLookup)]
    numFailedCodes <- length(failedCodes) 
    
    #count num successful matches
    numMatchedCountries <- nrow(dF) - numFailedCodes
    #printing info to console    
    cat(numMatchedCountries,"codes from your data successfully matched countries in the map\n")
           

    #failedCountries : reports on names of failed countries 
    #if user has specified the name of a country column in the function call
    failedCountries <- dF[[nameCountryColumn]][is.na(matchPosnsInLookup)]
    failedCountries <- cbind(failedCodes,"failedCountries"=as.character(failedCountries))
    
    #printing info to console    
    cat(numFailedCodes,"codes from your data failed to match with a country code in the map\n")
    if (verbose) print(failedCountries)
    
 #     failedCodes failedCountries                 
 #[1,] "CIV"       "Ivory Coast"                   
 #[2,] "COD"       "Congo, Democratic Republic"    
                       
    #put something here to try to match failed countries or codes 
    #initially just for if country name used as the join 
    #if ( suggestForFailedCodes )
    #if ( suggestForFailedCodes && joinCode=="NAME")
    #   {
    #    for( i in 1 : numFailedCodes)
    #       {
    #        correctCountry <- getCountryName(failedCodes[i])
    #        if (!is.na(correctCountry))
    #           {
    #            #dF[[joinCode]][is.na(matchPosnsInLookup)]
    #            dF[[joinCode]][which(dF[[joinCode]]==failedCodes[i])] <- correctCountry
    #           }          
    #       }
    #   }
    
    #can also get at countries in the lookup that don't appear in user data, by reversing match arguments
    matchPosnsInUserData <- match(as.character(mapWithData@data[[joinCode]])
                                , as.character(dF[[joinCode]])) 
    
    #these are the codes in lookup that aren't found in user data
    codesMissingFromUserData <- as.character( mapWithData@data[[joinCode]][is.na(matchPosnsInUserData)] )                            
    countriesMissingFromUserData <- as.character( mapWithData@data[["NAME"]][is.na(matchPosnsInUserData)] )
      
    numMissingCodes <- length(codesMissingFromUserData) 
    
    #printing info to console
    cat(numMissingCodes,"codes from the map weren't represented in your data\n")
    #if (verbose) #if (verbose) print more messages to console
    #   {
    #    if (nameJoinColumn!="NAME")                             
    #    {   print(cbind(codesMissingFromUserData,countriesMissingFromUserData))
    #    }else #if joined on NAME don't want to print names twice 
    #        print(codesMissingFromUserData)
    #   } # 


    #merging lookup table onto user data for those codes that match
    #dF2 <- cbind(dFlookupCodes[matchPosnsInLookup,],dF)    
    #the other way around to before, i.e. joining data onto map
    
    mapWithData@data <- cbind(mapWithData@data, dF[matchPosnsInUserData,])

    #test colouring map by region & subregion seems to show order has been retained
    #plot(mapWithData,col=mapWithData@data$REGION)
    #plot(mapWithData,col=mapWithData@data$SUBREGION)

    #returning the sPDF with the user data joined to the map polygons
    invisible(mapWithData)
   
   } #end of joinCountryData2Map()

#testing these should give perefct matches
#sPDF <- joinCountryData2Map(tlow,joinCode='NAME',nameJoinColumn='ADMIN',verbose=T,mapResolution="low")
#sPDF <- joinCountryData2Map(tlow,joinCode='NAME',nameJoinColumn='ADMIN',verbose=T,mapResolution="coarse")
#data(countryExData)
#sPDF <- joinCountryData2Map(countryExData,joinCode='NAME',nameJoinColumn='Country',verbose=T,mapResolution="low")

#[12,] "United States"                    "United States"                   
#[13,] "Viet Nam"                         "Viet Nam"                        
#106 codes from the map weren't represented in your data
#> rwmGetISO3("Viet Nam")
#[1] "VNM"

