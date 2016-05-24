### Parse Morningstar spreadsheets

# This is a parser for a spreadsheet distributed by Morningstar for their 
# Morningstar Direct product.  The spreadsheet includes a method for downloading
# historical data, whether daily, weekly, or monthly.  

# The function is designed to read in the spreadsheet and parse the columns into
# either separate columns, which it stores as instruments, or as a block, which
# it passes back to be assigned.  

# The "peers" method is for importing returns reported for hedge funds or other
# fund types, typically for use with PerformanceAnalytics.  

# The "instruments" method is for importing returns for a list of indexes, which
# are converted go CR objects ( Close, Return), stored in the specified directory, 
# and registered as instruments using FinancialInstrument.

## Load packages
parse.Morningstar <- function (filename, make.instruments=FALSE, filesroot = "~/Data/Morningstar", sheetname = "Calendar Return", symbol.row=c("Symbol", "SecID", "ISIN", "CUSIP", "Name")) {
  require(PerformanceAnalytics)
  require(gdata)
  require(quantmod)

  symbol.row = symbol.row[1]

  x = read.xls(filename, sheet=sheetname, check.names=FALSE, stringsAsFactors=FALSE, header=FALSE, strip.white=TRUE)
  ISOdates = as.Date(as.yearmon(x[-1:-7,1], format="%m/%Y"), frac=1)
  x.xts = as.xts(data.matrix(x[-1:-7,-1])/100, order.by=ISOdates)
  # mangles ampersands when read in, so fix them in colnames
  x.colnames=gsub("amp;", "", x[7, -1],fixed=TRUE)
  colnames(x.xts)=x.colnames
  
  if(!make.instruments) # just return the object to the work space
    return(x.xts)
  else{
    # Parse and create individual instruments with FI
    
    # Set the working directory to the folder that contains
    # the resulting data.
    setwd(filesroot)
    
    # Create symbol names
    # Identify which row to use for symbol name
    row = grep(symbol.row, x[3:7,1], value=FALSE) + 2
    symbolNames = x[row,-1]
    if (any(symbolNames=="")){
      warning(paste("The '", symbol.row, "' field has missing values.  Falling back to use 'SecId' instead.", sep=""), immediate.=TRUE)
      symbolNames = x[3,-1]
    }
    
    for(i in 1:length(symbolNames)) {
      # check to make sure directories exist for each symbol
      dir.create(paste(filesroot, symbolNames[i], sep="/"), showWarnings = FALSE, 
      recursive = FALSE, mode = "0777")
    }
    
    # Parse the columns into individual return objects
    print("Processing columns as symbols...")
    for( i in 1:dim(x.xts)[2]){
      # Column
      tmp = x.xts[,i]
      # Create an index from the returns
      tmp.desc = colnames(tmp)
      colnames(tmp)="Returns"
      xtsAttributes(tmp) <- list(Description = tmp.desc)
      save(tmp, file=paste(filesroot, symbolNames[i], paste(symbolNames[i], ".rda", sep=""), sep="/"))
      print(paste(symbolNames[i],", ", tmp.desc, sep=""))
      instrument(symbolNames[i], currency="USD", multiplier=1, tick_size=.0001, start_date=head(na.omit(index(tmp)),1), description=tmp.desc, data="R", source="Morningstar", assign_i=TRUE)
    }
    
    # Now, whenever you log in you need to register the instruments.  This
    # might be a line you put into .Rprofile so that it happens automatically:
    # require(quantmod) # this requires a development build after revision 560 or so.
    setSymbolLookup.FI(base_dir=filesroot, split_method='common')
    print( "Now, whenever you log in you need to register the instruments.  This")
    print( "might be a line you put into .Rprofile so that it happens automatically:")
    print( "> require(quantmod) # this requires a dev build after revision 560 or so.")
    print( paste("> setSymbolLookup.FI(base_dir=",filesroot,", split_method='common')"))
    print( "Now you should be able to use getSymbols to retrieve your data.")
  
  }
}
# Now you should be able to:
# getSymbols("FOUSA08MS7")
# or whatever...
# attr(FOUSA08MS7, "Description")
# 
# charts.PerformanceSummary(FOUSA08MS7, ylog=TRUE, wealth.index=TRUE, main = attr(FOUSA08MS7, "Description"))
# tail(FOUSA08MS7)

# symbolNames
# [1] "SecId"      "FOUSA08MS7" "FOUSA08MS6" "FOUSA08MS8"
# Unfortunately, these are the identifiers for the S&P GSCI Reduced Energy Index

# > attr(FOUSA08MS6, "Description")
# [1] "S&P GSCI Reduced Energy Spot"
# > attr(FOUSA08MS8, "Description")
# [1] "S&P GSCI Reduced Energy TR"
# > attr(FOUSA08MS7, "Description")
# [1] "S&P GSCI Reduced Energy Excess Return"
# > roll.R=FOUSA08MS7-FOUSA08MS6
