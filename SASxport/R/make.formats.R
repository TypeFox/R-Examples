## Take a list of data frames and process factor objects:
##
## For each factor object
##   1) generate SAS format information
##   2) add the "factor" attribute to these factors with the name of the generated SAS format
##   3) add this SAS format information to the FORMATS dataframe, creating it if necessary.
##
## Then return a new list of dataframe containing
##   1) The (potentially modified) data frames
##   2) the FORMAT dataframe
##
## Assumptions:
##   - dfList is a list of data frames
##   - an element named FORMAT format, if it exists, contains SAS format information
##   - each element of the list is named
##
make.formats <- function( dfList, formats=NULL )
  {
    if(missing(formats) || is.null(formats) )
      formats <- dfList$FORMATS
    
    dfList$FORMATS <- NULL
    
    if(is.null(formats))
      formats <- empty.format.table()

    retval <- list()
    formatIndex <- 1
    for(dfName in names(dfList))
      {
        # get the df we're working on
        df <- dfList[[dfName]]

        for(varName in colnames(df))
          {
            var <- df[[varName]]
            if(is.factor(var) && is.null(SASformat(var)) )
              {
                # We need unique format names, but SAS restricts
                # format names alpha characters.  To create a unique
                # name, convert counter to digits, then convert digits
                # (0-9) to letters (A-J)
                # Note: this mechanism will fail if more than 10,000
                #       factor names are needed...
                alphaStr <- chartr("0-9","A-J",as.character(formatIndex-1))
                formatName <-  paste("RFMT", alphaStr, sep="")
                formatIndex <- formatIndex+1
                
                formats <- rbind(formats,
                                 make.format.factor(var, formatName )
                                 )
                
                SASformat(var) <- formatName
                df[[varName]] <- var
            }
          }
        dfList[[dfName]] <- df
        
        # copy the df over to the return object
        retval[[dfName]] <- df

      }


    retval$FORMATS <- formats
    
    return( retval );

  }


make.format.factor <- function(var, fname)
  {
    if(missing(fname))
      formatName <- make.names(deparse(substitute(var)))
    else
      formatName <- fname
    varLevels <- levels(var)
    formats <- empty.format.table()

    if(nlevels(var)>0)
      for( j in 1:nlevels(var) )
        formats <- add.format.entry(formats,
                                    formatName,
                                    j,
                                    j,
                                    varLevels[j]
                                    )
    formats
  }


empty.format.table <- function()
  {
    formats <- data.frame(
                          FMTNAME = character(0), 
                          START = character(0), 
                          END = character(0), 
                          LABEL = character(0), 
                          MIN = integer(0), 
                          MAX = integer(0), 
                          DEFAULT = integer(0), 
                          LENGTH = integer(0), 
                          FUZZ = integer(0), 
                          PREFIX = character(0), 
                          MULT = integer(0), 
                          FILL = character(0), 
                          NOEDIT = integer(0), 
                          TYPE = character(0), 
                          SEXCL = character(0), 
                          EEXCL = character(0), 
                          HLO = character(0), 
                          DECSEP = character(0), 
                          DIG3SEP = character(0), 
                          DATATYPE = character(0), 
                          LANGUAGE = character(0)
                          )
  }

add.format.entry <- function(
                             formats,
                             FMTNAME, 
                             START,
                             END,
                             LABEL,
                             MIN = 1,
                             MAX = 40,
                             DEFAULT = 6,
                             LENGTH = 6, 
                             FUZZ = 1e-12,
                             PREFIX = "",
                             MULT = 0,
                             FILL = "",
                             NOEDIT = 0,
                             TYPE = "N",
                             SEXCL = "N", 
                             EEXCL = "N",
                             HLO = "",
                             DECSEP = "",
                             DIG3SEP = "",
                             DATATYPE = "",
                             LANGUAGE = ""
                             )
  {
    rbind(formats,
          data.frame(
                     FMTNAME = as.character(FMTNAME), 
                     START   = as.character(START), 
                     END     = as.character(END), 
                     LABEL   = as.character(LABEL), 
                     MIN     = as.integer(MIN), 
                     MAX     = as.integer(MAX), 
                     DEFAULT = as.integer(DEFAULT), 
                     LENGTH  = as.integer(LENGTH), 
                     FUZZ    = as.integer(FUZZ), 
                     PREFIX  = as.character(PREFIX), 
                     MULT    = as.integer(MULT), 
                     FILL    = as.character(FILL), 
                     NOEDIT  = as.integer(NOEDIT), 
                     TYPE    = as.character(TYPE), 
                     SEXCL   = as.character(SEXCL), 
                     EEXCL   = as.character(EEXCL), 
                     HLO     = as.character(HLO), 
                     DECSEP  = as.character(DECSEP), 
                     DIG3SEP = as.character(DIG3SEP), 
                     DATATYPE = as.character(DATATYPE), 
                     LANGUAGE = as.character(LANGUAGE)
                     )
          )
  }
