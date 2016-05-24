#' Internal function allowing map extents to be set from area names
#' 
#' Allows map extents to be set from country or area names (e.g. India, Africa
#' )
#' 
#' Can be called by \code{\link{mapCountryData}} and
#' \code{\link{mapGriddedData}}
#' 
#' @param mapRegion a country name from getMap()[['NAME']] or one of
#' 'eurasia','africa','latin america','uk','oceania','asia'
#' @return a dataframe containing we,ea,so,no values in degrees between -180 &
#' +180
#' @author andy south
#' @keywords dplot
#' @examples
#' 
#' mapCountryData( mapRegion='Africa' )
#' mapCountryData( mapRegion='India' )
#' 
#' @export setMapExtents
setMapExtents <- function(mapRegion='world')
   {  
    #sets map region from names passed to it
    #by returning a data frame containing wesn
    #maybe also define ocean regions
    
    
    #adding ability to set mapRegion from any countries
    map <- getMap()
    if ( mapRegion %in% map@data$ADMIN )
       {
        bb <- bbox(map[ map@data$ADMIN == mapRegion, ])
        dFmapExtents <- data.frame( we=bb[1], ea=bb[3], so=bb[2], no=bb[4] )
        return(dFmapExtents)       
       }
    
    
    listMapRegions=c('eurasia','africa','latin america','north america','uk','oceania','asia')
    if ( mapRegion == 'world' )
       {}else
    if ( mapRegion == 'eurasia' | mapRegion == 'Eurasia' )#1
       {we=-20;   ea=110;   so=20;   no=70} else 
    if ( mapRegion == 'africa' | mapRegion == 'Africa' )#2
       {we=-20;   ea=55;    so=-35;   no=38} else 
    if ( mapRegion == 'latin america' | mapRegion == 'Latin America' )#3
       {we=-118;   ea=-35;    so=-58;   no=30} else 
    if ( mapRegion == 'north america' | mapRegion == 'North America' )#3
       {we=-150;   ea=-70;    so=10;   no=85} else        
    if ( mapRegion == 'uk' | mapRegion == 'UK' )#4
       {we=-10;   ea=5;    so=50;   no=70} else                  
    if ( mapRegion == 'oceania' | mapRegion == 'Oceania' )#5
       {we=50;   ea=180;    so=-50;   no=0} else 
    if ( mapRegion == 'asia' | mapRegion == 'Asia' )#6
       {we=60;   ea=140;    so=-15;   no=55} else
    if ( mapRegion == 'europe' | mapRegion == 'Europe' )#1
       {we=-10;   ea=45;   so=35;   no=70} else        
       { 
        warning("The mapRegion you specified(",mapRegion,") needs to be a country name from getMap()$NAME or one of : (",paste(listMapRegions,""),").Plotting whole world") 
        we=-160;ea=160;so=-80;no=90
       }        

    dFmapExtents <- data.frame(we=we, ea=ea, so=so, no=no )
    
    return(dFmapExtents)
    
    } #end of setMapExtents function

