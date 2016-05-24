ConvertMedUnits <- function(x, measurement, abbreviation,
                            to=c("Conventional","SI","US"),
                            exact=!missing(abbreviation))
  {
    MedUnits <- NULL ## Define to avoid R CMD check warning
    data(MedUnits,package='gdata', envir=environment())
    
    to=match.arg(to)
    if(!missing(measurement) && missing(abbreviation))
      {
        if(exact)
          matchUnits <- MedUnits[tolower(MedUnits$Measurement)==
                                 tolower(measurement),]
        else
          matchUnits <- MedUnits[grep(measurement, MedUnits$Measurement,
                                  ignore.case=TRUE),]
      }
    else if(missing(measurement) && !missing(abbreviation))
      {
        if(exact)
          matchUnits <- MedUnits[tolower(MedUnits$Abbreviation)==
                                 tolower(abbreviation),]
    else
      matchUnits <- MedUnits[grep(match, MedUnits$Abbrevation,
                                  ignore.case=TRUE),]
      }
    else # both missing or both specified
      stop("One of `measurement' or `abbreviation' must be specified.")


    if(nrow(matchUnits)>1)
      stop(
           paste("More than one matching row.  Please use 'exact=TRUE' ",
                 "and supply one of these matching strings:",
                 paste('\t"',matchUnits$Measurement, '"', sep='', collapse="\n\t"),
                 sep="\n\t"))
   else if (nrow(matchUnits)<1)
     stop("No match")

    if (to %in% c("Convetional", "US"))
      {
        retval <- x / matchUnits$Conversion
        attr(retval,"units") <- matchUnits$ConventionalUnits
      }
    else
      {
        retval <- x * matchUnits$Conversion
        attr(retval,"units") <- matchUnits$SIUnits
      }
    retval
  }
      
    

    
    
