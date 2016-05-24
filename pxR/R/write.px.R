#################################################################
# 
# File:         write.px.R
# Purpose:      Write an object of class 'px' to a PC-Axis file 
#
# Created:      20110618
# Authors:      fvf, cjgb, opl
#
# Modifications: 
#               fvf (20130618)
#               cjgb, 20130811: added on.exit hook to close the open connection
#               cjgb, 20130811: fixed encoding issues (testing)
#               cjgb, 20130813: we do not generate files using KEYS (i.e., if present,
#                         the KEYS part needs to be ignored)
#               cjgb, 20131117: some users want to specify heading & stub precisely
#                         in order to control the shape of the resulting matrix
#################################################################


write.px <- function ( obj.px, filename, heading = NULL, stub = NULL, fileEncoding = "ISO-8859-1" )
{
  
  if ( ! inherits( obj.px, "px" ) )
    stop("Error: object needs to have class 'px'")
  
  ## auxiliary functions ## 
  unquote <- function(x)
    gsub( '^"|"$', "", x)           # removes " at beginning/end of string
  
  requote <- function(x)            # adds dquotes (") to character vectors
    paste( '"', unquote(x), '"', sep = "")
    
  wf  <- function ( ..., my.encoding = ifelse(fileEncoding == "ISO-8859-1", "latin1", "utf8") ) {
    cadena <- paste(..., sep = "")
    cadena <- iconv(cadena, my.encoding)
    cat( cadena, file = con, sep = "" ) 
  }
  ## end - auxiliary functions ##
  
  # modify px object providing sensible default values
  obj.px[['LAST.UPDATED']]$value <- format(Sys.time(), "%Y%m%d %H:%M:%S")
  obj.px[['CHARSET']]$value      <- ifelse( fileEncoding == "ISO-8859-1", 'ANSI', fileEncoding )
  obj.px$INFO$value              <- "File generated using R and package pxR (http://pxr.r-forge.r-project.org/)"
  
  # we do not generate files that use the KEYS argument
  obj.px$KEYS <- NULL
  
  # obj.px names may have changed (and got - changed into . by R)
  # we need to revert that  
  names(obj.px) <- gsub("\\.", "-", names(obj.px))
  
  # we want to write the output file keywords in the 'right' order
 
  order.kw <- c("CHARSET", "AXIS-VERSION", "CODEPAGE", "LANGUAGE",
                "LANGUAGES", "CREATION-DATE", "NEXT-UPDATE", "PX-SERVER",
                "DIRECTORY-PATH", "UPDATE-FREQUENCY", "TABLEID", "SYNONYMS",
                "DEFAULT-GRAPH", "DECIMALS", "SHOWDECIMALS", "ROUNDING",
                "MATRIX", "AGGREGALLOWED", "AUTOPEN", "SUBJECT-CODE",
                "SUBJECT-AREA", "CONFIDENTIAL", "COPYRIGHT", "DESCRIPTION",
                "TITLE", "DESCRIPTIONDEFAULT", "CONTENTS", "UNITS", "STUB",
                "HEADING", "CONTVARIABLE", "VALUES", "TIMEVAL", "CODES",
                "DOUBLECOLUMN", "PRESTEXT", "DOMAIN", "VARIABLE-TYPE",
                "HIERARCHIES", "HIERARCHYLEVELS", "HIERARCHYLEVELSOPEN",
                "HIERARCHYNAMES", "MAP", "PARTITIONED", "ELIMINATION", "PR",
                "ECISION", "LAST-UPDATED", "STOCKFA", "CFPRICES", "DAYADJ",
                "SEASADJ",  "CONTACT", "REFPERIOD", "BASEPERIOD",
                "DATABASE", "SOURCE", "SURVEY", "LINK", "INFOFILE",
                "FIRST-PUBLISHED", "META-ID", "OFFICIAL-STATISTICS", "INFO",
                "NOTEX", "NOTE", "VALUENOTEX", "VALUENOTE", "CELLNOTEX",
                "CELLNOTE", "DATASYMBOL1", "DATASYMBOL2", "DATASYM", "BOL3",
                "DATASYMBOL4", "DATASYMBOL5", "DATASYMBOL6", "DATASYMBOLSUM",
                "DATASYMBOLNIL", "DATANOTECELL", "DATANOTESUM", "DATANOTE",
                "KEYS", "ATTRIBUTE-ID", "ATTRIBUTE-TEXT", "ATTRIBUTES",
                "PRECISION","DATA")
 
  order.px  <- charmatch(names(obj.px), order.kw, nomatch=999)   # keyword no in list to end
  new.order <- setdiff( names(obj.px)[order(order.px)], "DATA" ) # all but "DATA"
  
  if(! is.null(heading)){
    if(is.null(stub))
      stop("If heading is specified, you need also specify the stub parameter.")
    
    if(! setequal(c(heading, stub), c(obj.px$HEADING$value, obj.px$STUB$value)) )
      stop("Specified heading and stub parameters differ from those in the px object")
    
    obj.px$HEADING$value <- heading
    obj.px$STUB$value    <- stub
  }
 
  ## open the connection and close automatically on exit
  con <- file( description = filename, open = "w", encoding = fileEncoding )
  on.exit(close(con))
  
  ## write new file 
  #cat('', file = con, append = F) # init file
  
  ## metadata part
  for (key in new.order ) {
    
    # these are exceptions: no quotes
    # e.g.: 'DECIMALS=0;'
    if( key %in% c('DECIMALS', 'SHOWDECIMALS', 'ELIMINATION', 'COPYRIGHT', 'DESCRIPTIONDEFAULT', 'DAYADJ', 'SEASADJ')){
      wf( key, "=")
      wf( unquote(obj.px[[key]]$value) )
      wf(';\n')
      next
    }
    
    # most metadata entries are like this: 
    # 'KEY="a","b";'
    if ( names(obj.px[[key]])[1] == 'value' ) {  
      wf( key, "=")
      wf( paste( requote(obj.px[[key]]$value), collapse = ',') )
      wf(';\n')
      next
    } 
    
    # meta with second name; there can be more than one
    for (subkey in names(obj.px[[key]])){
      wf ( key, '("', subkey, '")=' ) 
      wf ( paste( requote( obj.px[[key]][[subkey]] ), collapse = ',') )
      wf ( ';\n' )              
    }     
  }
 
  # DATA part
  wf('DATA=\n')
  
  zz <- formatC(as.array(obj.px),
                format = 'f',
                digits = as.numeric(obj.px$DECIMALS$value),
                drop0trailing = T, flag = '-')
  #zz <- str_trim(zz) 
  zz <- gsub("NA", '".."', zz)
  
  zz <- aperm(zz, c(rev(obj.px$HEADING$value), rev(obj.px$STUB$value)))
  write(zz, file = con, ncolumns = sum( dim(zz)[1:length(obj.px$HEADING$value)] ), append = T )
  
  wf(";\n")
}  

