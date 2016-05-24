#' Get a table of trait measurements from YouTheria
#' 
#' Retrieves a \code{data.frame} of trait measurements with facilities to select by
#' location, species name and/or measurement type.
#' 
#' @param measurementType Measurement types to collect data for. If \code{NULL} (default),
#'        all measurement types are returned. Can also be 'numeric' or 'character' (or a
#'        list of either type) and will filter by MeasurementTypeID and the measurement
#'        respectivly. MeasurementTypeIDs and names can be found using getMeasurementTypes().
#' @param MSW93Binomial Character giving the latin name of a species (or list of species)
#'        for which measurements are required. Naming should follow Mammal Species of the 
#'        World 1993.
#' @param MSW05Binomial Character giving the latin name of a species (or list of species)
#'        for which measurements are required. Naming should follow Mammal Species of the 
#'        World 2005.
#' @param country Character specifying the country from which you wish to collect data.
#'        If \code{NULL} all data is retrieved. If specified then \code{locationOnly}
#'        is set to \code{TRUE}
#' @param StudyUnitId Numeric specifying the StudyUnitId from which you wish to collect data.
#'        If \code{NULL} all data is retrieved
#' @param locationData Logial dictating whether location information should be added to the
#'        output. Defualt is \code{FALSE} but is set to \code{TRUE} if either StudyUnitId or
#'        country are specified.
#' @param locationOnly If \code{TRUE} data is only be returned if it has location information.
#' @param cast If \code{TRUE} (default) then the data is cast so that each observation 
#'        is one row in the output. If false then each observation has one row for 
#'        each data element recorded (i.e. range, mean, units, etc)
#' @param silent If \code{TRUE} progress reporting is silenced
#' 
#' @return A \code{data.frame} with each row giving a trait measurement        
#' @export
#' @import reshape2
#' @examples
#' 
#' \dontrun{
#' # Select measurement type by id
#' M14 <- getMeasurementData(14)
#' M22_7_2 <- getMeasurementData(c(22,7,2))
#' 
#' # Select measurement type by name
#' WM <- getMeasurementData('Wing Morphology')
#' WM_TN <- getMeasurementData(c('Wing Morphology','Teat Number'))
#' 
#' # Select by measurement type and species name
#' PpPr_bodymass <- getMeasurementData(measurementType = 1,
#'                                     MSW93Binomial = c('Pongo pygmaeus','Peroryctes raffrayana'))
#'                    
#' #Select by measurement type, species name and location
#' Ob_Activity_Tanz <- getMeasurementData(measurementType = 'Activity Cycle',
#'                                        MSW05Binomial = 'Oryx beisa',
#'                                        country = 'Tanzania')
#' }

getMeasurementData <-
  function(measurementType = NA,
           MSW93Binomial = NA,
           MSW05Binomial = NA,           
           country=NULL,
           StudyUnitId=NULL,
           locationData = TRUE,
           locationOnly = FALSE,
           cast = TRUE,
           silent = FALSE
           ){
    # Prevent download of everthing (this will most likley cause a crash)
    if(all(is.na(c(measurementType, MSW93Binomial, MSW05Binomial))) &
       all(is.null(c(country, StudyUnitId)))) {
      stop('Downloading everything will crash your computer, use one of the filtering arguments')
    }
    
    if(!identical(measurementType,NA)){
            
        if(is.numeric(measurementType)){
          if(!exists('MTs')) MTs <- getMeasurementTypes()
          bad <- measurementType[!measurementType %in% MTs$Id]
          if(length(bad)==length(measurementType)){
            stop('All measurement types specified are invalid: Measurement Type ID(s) ',bad,' is/are not known. Use getMeasurementTypes() to find out what is appropriate.', sep='')
          } else if(length(bad)!=0){
            warning('Some measurement types unknown: Measurement Type ID(s) ',bad,' is/are not known. Use getMeasurementTypes() to find out what is appropriate.', sep='')
          }         
        }
        
        if(is.character(measurementType)){
          if(!exists('MTs')) MTs <- getMeasurementTypes()
          bad <- measurementType[!measurementType %in% MTs$Name]
          if(length(bad)==length(measurementType)){
            stop('All measurement types specified are invalid: Measurement Type ID(s) ',bad,' is/are not known. Use getMeasurementTypes() to find out what is appropriate.', sep='')
          } else if(length(bad)!=0){
            warning('Some measurement types unknown: Measurement Type ID(s) ',bad,' is/are not known. Use getMeasurementTypes() to find out what is appropriate.', sep='')
          }           
          measurementType <- MTs$Id[MTs$Name %in% measurementType]
        }
    }
    
    if(!identical(MSW93Binomial,NA) & !identical(MSW05Binomial,NA)) stop ('Cannot filter by MSW05Binomial and MSW93Binomial. Choose one or the other')
    
    if(!is.null(StudyUnitId)&!is.null(country)) stop('Cannot use both StudyUnitId and country at the same time')
   
    # This rather complicated bit of code creates all possible combinations
    # for the URL string
    URL <- apply(expand.grid(paste('?id=',
                                   apply(expand.grid(as.character(ifelse(is.na(measurementType),'',measurementType)),
                                                     ifelse(is.na(MSW93Binomial),'',MSW93Binomial)),
                                         1,paste,collapse='&MSW93Binomial='),
                             sep=''),
                             ifelse(is.na(MSW05Binomial),'',MSW05Binomial)),
                 1,paste,collapse='&MSW05Binomial=')

    if(!is.null(country) | !is.null(StudyUnitId)) locationData = TRUE
    
    if(locationData){
      
      if(!silent) cat('Retrieving trait information...')
      out <- runURL(URL,'m')
      if(!silent) cat('complete\n')
      
      # cast the data if requested
      if(length(out) == 0){        
        warning('No data was returned. Check species names are correct')
        out <- NULL        
      } else {
        
        if(cast & nrow(out) > 0){
          if(!silent) cat('Casting data...')
          out <- dcast(out,  MeasurementTypeID+MeasurementSetID+StudyUnitId+Genus+Species+SubSpecies
                       +MSW93Binomial+MSW05Binomial+AuthorityText ~ ValueType, value.var = 'MValue')    
          if(!silent) cat('complete\n')
        }
        
        if(!silent) cat('Retrieving location information...')
        loc_data <- getLocData(country, StudyUnitId)
        if(!silent) cat('complete\n')
        
        if(!silent) cat('Combining data...')
        startRowCount <- nrow(out)
        if(is.null(country) & is.null(StudyUnitId)){
          out <- merge(out, loc_data, all.x = !locationOnly, all.y = F)
        } else {
          out <- merge(out, loc_data, all.x = FALSE, all.y = F)
        }
        endRowCount <- nrow(out)
        if(!silent) cat('complete\n')
        
        if(locationOnly & startRowCount != endRowCount){
          if(endRowCount == 0){
            stop('No measurements have location information')
          } else {
          warning(paste(startRowCount - endRowCount, 'rows have been removed as they are missing location information'))
          }
        }
      } 
    } else {
      if(!silent) cat('Retrieving trait information...')
      out <- runURL(URL,'m')
      if(!silent) cat('complete\n')
      
      # cast the data if requested
      if(length(out) == 0){
        warning('No data was returned. Check species names and locations are correct')
        out <- NULL
      } else if(cast & nrow(out) > 0 & length(out) != 0){
        if(!silent) cat('Casting data...')
        out <- dcast(out,  MeasurementTypeID+MeasurementSetID+StudyUnitId+Genus+Species+SubSpecies
                     +MSW93Binomial+MSW05Binomial+AuthorityText ~ ValueType, value.var = 'MValue')    
        if(!silent) cat('complete\n')
      }
    }
    
    #Report species that are missing from the results be requested
    if(!identical(MSW93Binomial,NA)){
      missing <- MSW93Binomial[!MSW93Binomial %in% unique(out$MSW93Binomial)]
      if(length(missing)>0) warning(paste('There were no results returned for the following species:', paste(missing, collapse=', ')))    
    }
    if(!identical(MSW05Binomial,NA)){
      missing <- MSW05Binomial[!MSW05Binomial %in% unique(out$MSW05Binomial)]
      if(length(missing)>0) warning(paste('There were no results returned for the following species:', paste(missing, collapse=', ')))    
    }
    
    return(out)
  }
