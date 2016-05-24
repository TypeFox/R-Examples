download.DJUBS <- function (filesroot = "~/Data/DJUBS", download=TRUE) {
  # Script for parsing DJUBS index daily price data series from the
  # DJ website.
  
  # Peter Carl
  
  # DETAILS
  # Parse index close prices from the spreadsheet containing the full series:
  # http://www.djindexes.com/mdsidx/downloads/xlspages/ubsci/DJUBS_full_hist.xls
  
  
  # Several series, all index values
  # Remove the footer at the bottom
  # Load needed packages:
  require(xts)
  require(PerformanceAnalytics)
  require(gdata)
  require(FinancialInstrument)
  require(quantmod)
  currency("USD")
  # filesroot: Set the working directory, where there's a .incoming folder that 
  # contains the downloaded spreadsheet.
  
  # Create and set the working directory if it doesn't exist
  if (!file.exists(filesroot))
    dir.create(filesroot, mode="0777")
  
  # Create and set the .incoming directory if it doesn't exist
  if (!file.exists(paste(filesroot, "/.incoming", sep="")))
    dir.create(paste(filesroot, "/.incoming", sep=""), mode="0777")
  setwd(paste(filesroot, "/.incoming", sep=""))
  
  # Remove the old file from .incoming
  if(file.exists("DJUBS_full_hist.xls"))
    if(download){
      system("rm DJUBS_full_hist.xls")
  
      # Download the xls workbook directly from the web site:
      print("Downloading excel spreadsheet from DJUBS web site...")
      system("wget http://www.djindexes.com/mdsidx/downloads/xlspages/ubsci_public/DJUBS_full_hist.xls")
    }
  
  if(!file.exists("DJUBS_full_hist.xls"))
    stop(paste("No spreadsheet exists.  Download the spreadsheet to be processed from www.djindexes.com into ", filesroot, "/.incoming", sep=""))
  
  sheetnames=c("Excess Return", "Total Return")
  for(sheet in sheetnames){  
    print(paste("Reading", sheet, "sheet... This will take a moment..."))
    x = read.xls("DJUBS_full_hist.xls", sheet=sheet, stringsAsFactors=FALSE)
    
    # Add column names, get the descriptions to add as attributes
    colnames(x)=t(as.data.frame(apply(x[2,], FUN=as.character, MARGIN=1), stringsAsFactors=FALSE))
    x.attr = t(as.data.frame(x[1,], stringsAsFactors=FALSE))
    x=x[-1:-2,]
    
    # Get rid of the last line, which contains the disclaimer
    x=x[-dim(x)[1],]
    # Remove blank columns between sections
    x=x[,-which(apply(x,2,function(x)all(is.na(x))))]
    
    # Get attributes and labels
    categoryNames = x.attr[!is.na(x.attr)]
    symbolNames = paste(make.names(colnames(x[,])), ".IDX", sep="")
    symbolNamesMonthly = paste(make.names(colnames(x[,])), ".M.IDX", sep="")
    ISOdates = as.Date(x[,1], "%m/%d/%Y")
    
    for(i in 2:length(symbolNames)) {
      # check to make sure directories exist for each symbol, first for daily series...
      dir.create(paste(filesroot, symbolNames[i], sep="/"), showWarnings = FALSE, 
      recursive = FALSE, mode = "0777")
      # ... then for monthly series
      dir.create(paste(filesroot, symbolNamesMonthly[i], sep="/"), showWarnings = FALSE, 
      recursive = FALSE, mode = "0777")
    }
    
    # Parse the columns into individual price objects
    print("Processing columns as symbols...")
    for( i in 2:dim(x)[2]){
      x.xts = as.xts(as.numeric(x[,i]), order.by=ISOdates)
      R.xts = Return.calculate(x.xts)
      x.xts = cbind(x.xts, R.xts)
      colnames(x.xts)=c(paste(symbolNames[i],".Close",sep=""), paste(symbolNames[i],".Returns",sep=""))
      xtsAttributes(x.xts) <- list(Description = paste(categoryNames[i], sheet, "Index"))
  
      save(x.xts, file=paste(filesroot, symbolNames[i], paste(symbolNames[i], ".rda", sep=""), sep="/"))
      print(paste(symbolNames[i],", ",categoryNames[i], ", ", sheet, sep=""))
      
      # Describe the metadata for each index
      instrument(symbolNames[i], currency="USD", multiplier=1, tick_size=.01, start_date=head(index(x.xts),1), description=paste(categoryNames[i], "Index"), data="CR", source="DJUBS", frequency="Daily", assign_i=TRUE)
      
      # Construct a monthly series from the daily series
      x.m.xts = to.monthly(Cl(x.xts))
      x.m.xts = cbind(x.m.xts[,4], Return.calculate(x.m.xts[,4]))
      colnames(x.m.xts)=c(paste(symbolNames[i],".Close",sep=""), paste(symbolNames[i],".Returns",sep=""))
      # @ TODO Want to delete the last line off ONLY IF the month is incomplete
      if(tail(index(x.xts),1) != as.Date(as.yearmon(tail(index(x.xts),1)), frac=1)) {
        # That test isn't quite right, but its close.  It won't work on the first
        # day of a new month when the last business day wasn't the last day of 
        # the month.  It will work for the second day.
        x.m.xts = x.m.xts[-dim(x.m.xts)[1],]
      }
        
      # Index is set to last trading day of the month.  
      # Reset index to last day of the month to make alignment easier with other monthly series.  
      index(x.m.xts)=as.Date(index(x.m.xts), frac=1)
        
      xtsAttributes(x.m.xts) <- list(Description = paste(categoryNames[i], sheet, "Index"))
  
      save(x.m.xts, file=paste(filesroot, symbolNamesMonthly[i], paste(symbolNamesMonthly[i], ".rda", sep=""), sep="/"))
      print(paste(symbolNamesMonthly[i],", ",categoryNames[i], ", ", sheet, sep=""))
      # Describe the metadata for each index
      instrument(symbolNamesMonthly[i], currency="USD", multiplier=1, tick_size=.01, start_date=head(index(x.xts),1), description=paste(categoryNames[i], "Index"), data="CR", source="DJUBS", frequency="Monthly", assign_i=TRUE)
    }
  }
  
  # Now, whenever you log in you need to register the instruments.  This
  # might be a line you put into .Rprofile so that it happens automatically:
  # require(quantmod) # this requires a development build after revision 560 or so.
  setSymbolLookup.FI(base_dir=filesroot, split_method='common')
  print( "Now, whenever you log in you need to register the instruments.  This")
  print( "might be a line you put into .Rprofile so that it happens automatically:")
  print( "> require(quantmod) # this requires a dev build after revision 560 or so.")
  print( "> setSymbolLookup.FI(base_dir=filesroot, split_method='common')")
  print( "Now you should be able to type:")
  print( "> getSymbols('DJUBSTR.IDX') ")
}

# # Now you should be able to:
# getSymbols("DJUBSTR.IDX")
# # chartSeries(Cl(DJUBSTR.IDX), theme="white")
# charts.PerformanceSummary(DJUBSTR.IDX[,"Returns"], ylog=TRUE, wealth.index=TRUE, main = "DJUBS Total Returns Index Returns")
# tail(DJUBSTR.IDX)
# 
# 
# symbols=c('DJUBSTR.M.IDX','DJUBS.M.IDX','DJUBSSP.M.IDX')
# getSymbols(symbols)
# 
# # Look at how Garett does this
# x=cbind(DJUBSTR.M.IDX[,2],DJUBS.M.IDX[,2],DJUBSSP.M.IDX[,2])
# y=NULL; for(i in symbols) y=c(y,attr(get(i),"Description"))
# y[3]="DJUBS Spot Index"
# colnames(x)=y
# 
# # Get an inflation series from FRED
# getSymbols("CPIAUCSL",src="FRED") #load CPI for inflation
# inflation = ROC(CPIAUCSL,12,type="discrete")/12
# index(inflation) = as.Date(as.yearmon(index(inflation,1)), frac=1)
# realspot=DJUBSSP.M.IDX[,2]-inflation
# 
# charts.RollingPerformance(realspot, width=12)
# # Not much of a difference during this period of low inflation
# 
# 
# 
# # To calculate roll returns (shown on daily)
# roll.R = DJUBS[,"Returns"]-DJUBSSP[,"Returns"]
# chart.CumReturns(roll.R)