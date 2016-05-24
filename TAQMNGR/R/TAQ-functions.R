################################################################################
##
## File:    TAQ-functions.R
##
## Purpose: Functions for package 'TAQ'.
##
## Created: 2012.02.13
##
## Version: 2013.11.22
##
## Todo list:
##   
################################################################################

################################################################################
## Part 1) General functions
##  .package.name()                  Name of the package.
################################################################################

.package.name <- 
function()
{
  ##############################################################################
  ## Description:
  ##  Set the name of the package.
  ##
  ## Arguments:
  ##  NONE
  ##
  ## Value:
  ##  (character[1]) package name.
  ##############################################################################

  ## FUNCTION:
  
  #### Answer
  "TAQMNGR"
}
# ----------------------------------------------------------------------


################################################################################
## Part 2) MAIN functions
##  TAQ.CleanTickByTick()                  Clean tick by tick data.
################################################################################


#TAQ.ControlFiles <- 
#function(dirInput, dirOutput) 
#{
#  ##############################################################################
#  ## Description:
#  ##
#  ## Arguments:
#  ##  NONE
#  ##
#  ## Value:
#  ##  (character[2]) labels.
#  ##############################################################################
#
#  ## FUNCTION:  
#
#  #### Settings
#  dirInput  <- path.expand(dirInput[1])
#  dirOutput <- path.expand(dirOutput[1])
#  
#  #### dirInput
#  x1 <- dirInput
#  if ( !file.exists(x1) )
#  {
#    stop("Directory dirInput '", x1, "' does not exist")
#  }
#
#  #### dirOutput
#  x1 <- dirOutput
#  if ( !file.exists(x1) )
#  {
#    ind <- dir.create(path = x1, recursive = FALSE)
#    if (!ind)
#    {
#      stop("Unable to create dirOutput '", x1, "'")
#    } 
#  }
#  
#  #### Answer
#  .Call( 
#    name = "ControlFiles",
#    DirIn = as.character(dirInput[1]),
#    DirOut = as.character(dirOutput[1]), 
#    PACKAGE = .package.name() )
#}
# ------------------------------------------------------------------------------


TAQ.CleanTickByTick <- 
function(dirInput, dirOutput, 
  window = 80, deltaTrimmed = 0.10, granularity = 0.04, useCleaned = TRUE)
{
  ##############################################################################
  ## Description:
  ##  Clean tick by tick data.
  ##
  ## Arguments:
  ##  dirInput: (character[1]) name of the input directory.
  ##  dirOutput: (character[1]) name of the output directory.
  ##  window: (numeric[1]) window size (an integer > 3). 
  ##  deltaTrimmed: (numeric[1]) trimming proportion (within (0,1)). 
  ##  granularity: (numeric[1]) a positive parameter that avoids a zero variance 
  ##   in trimming.
  ##  useCleaned: (logical[1]) use previously cleaned data (if any)? 
  ##
  ## Value:
  ##  NONE
  ##############################################################################

  ## FUNCTION:  

  #### Settings
  dirInput     <- .sanitize.dirname(x = path.expand(dirInput[1]))
  dirOutput    <- .sanitize.dirname(x = path.expand(dirOutput[1]))
  window       <- window[1]
  deltaTrimmed <- deltaTrimmed[1]
  granularity  <- granularity[1]
  useCleaned   <- as.logical(useCleaned[1])
  
  #### dirInput
  x1 <- dirInput
  if ( !file.exists(x1) )
  {
    stop("FATAL ERROR: folder ", x1, " not found.")
  }

  #### dirOutput
  x1 <- dirOutput
  if ( !file.exists(x1) )
  {
    ind <- dir.create(path = x1, recursive = FALSE)
    if (!ind)
    {
      stop("FATAL ERROR: unable to create ", x1, ".")
    } 
  }
  
  #### window
  x1 <- window
  if ( x1 < 3 )
  {
    stop("FATAL ERROR: argument 'window' must be >= 3")
  }
  
  #### granularity
  x1 <- granularity
  if ( x1 < 0 )
  {
    stop("FATAL ERROR: argument 'granularity' must be > 0")
  }

  #### deltaTrimmed
  x1 <- deltaTrimmed
  if ( x1 <= 0 || x1 >= 1)
  {
    stop("FATAL ERROR: argument 'deltaTrimmed' must be within (0,1)")
  }
  
  #### Answer
  .Call( "CleanTickByTick",
    DirIn    = as.character(dirInput),
    DirOut   = as.character(dirOutput),
    Win      = as.integer(window),
    DeltaTr  = as.double(deltaTrimmed),
    Gran     = as.double(granularity),
    UseClean = as.integer(useCleaned),
    PACKAGE  = .package.name() )
}
# ------------------------------------------------------------------------------


TAQ.Report <- 
function(dirInput, symbol)
{
  ##############################################################################
  ## Description:
  ##  Report cleaned data for one symbol.
  ##
  ## Arguments:
  ##  dirInput: (character[1]) name of the directory including tick-by-tick data
  ##   previously cleaned.
  ##  symbol: (character[1]) the symbol reported.
  ##
  ## Value:
  ##  NONE
  ##############################################################################

  ## FUNCTION:  

  #### Settings
  dirInput <- .sanitize.dirname(x = path.expand(dirInput[1]))
  symbol   <- symbol[1]
  
  #### Argument check
  ## dirInput
  x1 <- dirInput
  if ( !file.exists(x1) )
  {
    stop("FATAL ERROR: folder ", x1, " not found")
  }
  
  #### Answer
  .Call("CleaningReport",
    DirIn         = as.character(dirInput),
    Symbol        = as.character(symbol),
    PACKAGE       = .package.name() )
}
# ------------------------------------------------------------------------------


TAQ.Aggregate <- 
function(dirInput, symbol, bin, useAggregated = TRUE)
{
  ##############################################################################
  ## Description:
  ##  Aggregate tick-by-tick data previously cleaned.
  ##
  ## Arguments:
  ##  dirInput: (character[1]) name of the directory including tick-by-tick data
  ##   previously cleaned.
  ##  symbol: (character[1]) .
  ##  bin: (numeric[1]) the bin size in seconds.
  ##  useAggregated: (logical[1]) use previously aggregated data (if any)? 
  ##
  ## Value:
  ##  NONE
  ##############################################################################

  ## FUNCTION:  

  #### Settings
  dirInput      <- .sanitize.dirname(x = path.expand(dirInput[1]))  
  bin           <- round(as.numeric(bin[1]))
  useAggregated <- as.logical(useAggregated[1])
  nSym          <- NROW(symbol)
  
  #### Argument check
  ## dirInput
  x1 <- dirInput
  if ( !file.exists(x1) )
  {
    stop("FATAL ERROR: folder ", x1, " not found")
  }
  ## symbol
  if ( nSym == 0 )
  {
    stop("FATAL ERROR: argument 'symbol' cannot be empty")
  }
  ## bin
  if ( bin <= 0 )
  {
    stop("FATAL ERROR: argument 'bin' must be positive")
  }
  
  #### Answer
  for (i in 1 : nSym)
  {
    .Call("Aggregate",
      DirIn         = as.character(dirInput),
      symbol        = as.character(symbol[i]),
      bin           = as.character(bin),
      useAggregated = as.integer(useAggregated),
      PACKAGE       = .package.name() )
  }
}
## ------------------------------------------------------------------------------


TAQ.Read <- 
function(dirInput, symbol, import = NULL, startDate, endDate, bin)
{
  ##############################################################################
  ## Description:
  ##  Aggregate tick-by-tick data previously cleaned.
  ##
  ## Arguments:
  ##  dirInput: (character[1]) name of the directory including tick-by-tick data
  ##   previously cleaned.
  ##  symbol: (character[1]) .
  ##  import: (character) Fields imported from the following list:
  ##   "FIRST": 
	##   "MIN":
	##   "MAX": 
	##   "LAST":
	##   "SIZE":
	##   "#TRADES":
	##   "VWAP":
  ##  startDate: (numeric[1]) 
  ##  endDate: (numeric[1]) 
  ##  bin: (numeric[1]) the bin size in seconds.
  ##
  ## Value:
  ##  (data.frame)
  ##############################################################################

  ## FUNCTION:  

  #### Settings
  dirInput     <- .sanitize.dirname(x = path.expand(dirInput[1]))
  symbol       <- symbol[1]
  startDate    <- as.character(startDate[1])
  endDate      <- as.character(endDate[1])
  bin          <- round(as.numeric(bin[1]))
  import       <- as.character(import)
  importDef    <- c("FIRST", "MIN", "MAX", "LAST", "SIZE", "#TRADES", "VWAP")
  
  #### Import
  if ( NROW(import) == 0 )
  {
    import <- importDef
  }
  else
  {
    x1 <- substr(x = importDef, start = 1, stop = 7)
    x2 <- substr(x = import,    start = 1, stop = 7)
    ind <- x1 %in% x2
    import <- importDef[ind]
  }
  if ( NROW(import) == 0 )
  {
    stop("FATAL ERROR: argument 'import' must include at least one of the following elements ", 
      importDef)
  }
    
  #### dirInput
  x1 <- dirInput
  if ( !file.exists(x1) )
  {
    stop("FATAL ERROR: folder ", x1, " not found")
  }
  
  #### startDate
  x1  <- startDate
  nch <- nchar(x1)
  if ( nch < 4 )
  {
    stop("FATAL ERROR: argument 'startDate' must be in the format yyyymmdd")
  }
  yearStart  <- substr(x = x1, start = 1, stop = 4)
  monthStart <- ifelse(nch >= 6, substr(x = x1, start = 5, stop = 6), "01")
  dayStart   <- ifelse(nch >= 8, substr(x = x1, start = 7, stop = 8), "01")
  
  #### endDate
  x1 <- endDate
  if ( nchar(x1) < 6 )
  {
    stop("FATAL ERROR: argument 'endDate' must be in the format yyyymmdd")
  }
  yearEnd  <- substr(x = x1, start = 1, stop = 4)
  monthEnd <- ifelse(nch >= 6, substr(x = x1, start = 5, stop = 6), 12)
  dayEnd   <- ifelse(nch >= 8, substr(x = x1, start = 7, stop = 8), 31)

  #### 
  if ( bin <= 0 )
  {
    stop("FATAL ERROR: argument 'bin' must be positive")
  }
  
  #### Call
  x <- .Call("Import",
    DirIn      = as.character(dirInput),
    symbol     = as.character(symbol),
    import     = as.character(import), 
    yearStart  = as.integer(yearStart),
    monthStart = as.integer(monthStart),
    dayStart   = as.integer(dayStart),
    yearEnd    = as.integer(yearEnd),
    monthEnd   = as.integer(monthEnd),
    dayEnd     = as.integer(dayEnd),
    bin        = as.integer(bin),
    PACKAGE    = .package.name() )
    
  #### Extract
  # x1 <- as.integer( paste(yearStart, monthStart, dayStart, sep = "") )
  # x2 <- as.integer( paste(yearEnd, monthEnd, dayEnd, sep = "") )
  # ind <- x1 <= x[, "DATE"] & x[, "DATE"] <= x2
  # x[ind, , drop = FALSE]
}
# ------------------------------------------------------------------------------


#################################################################################
## AUXILIARY FUNCTIONS
#################################################################################

.sanitize.dirname <- 
function(x) 
{
  ##############################################################################
  ## Description:
  ##  Sanitize dirname.
  ##
  ## Arguments:
  ##  (character) directory names.
  ##
  ## Value:
  ##  (character) directory names adjusted.
  ##############################################################################

  ## FUNCTION:

  #### Settings
  platform <- .Platform$OS.type
  file.sep <- "/"
  
  x <- gsub(x = x, pattern = "\\", replacement = file.sep, fixed = TRUE)

  #### Remove repeated .Platform$file.sep (if any)
  x <- gsub(x = x, pattern = paste("[", file.sep, "]+", sep = ""), 
    replacement = file.sep)
  
  #### Remove final .Platform$file.sep (if any)
  nch <- nchar(x)
  x1 <- substr(x = x, start = nch, stop = nch)
  ind <- x1 == file.sep
  x[ind] <- substr(x = x[ind], start = 1, stop = nch - 1)
  
  if (platform == "windows")
  {
  	x <- gsub(x = x, pattern = file.sep, replacement = "\\", fixed = TRUE)
  }
  x
}
# ------------------------------------------------------------------------------
