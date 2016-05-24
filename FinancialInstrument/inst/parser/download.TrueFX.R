#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
## cores should not be greater than the number of gigs of available memory.
cores <- if(length(args) > 0L) { as.numeric(args[1]) } else 4L

# Download Bid/Ask tick data for 15 FX pairs since May 2009 from TrueFX.com
#
# Garrett See 
#
# You probably _should_ sign up for a FREE account with TrueFX before 
# downloading their data, although presently (2012-04), it works without an 
# account.
#
# Data will be downloaded, unzipped, converted to xts, and saved as both
# tick data and 1 second frequency data.  Then it is saved by day.  
# seealso ?saveSymbols.days
#
# 3 directories will be created underneath base_dir: "archive", "tick", "sec"
#  - The archive directory will contain the original zip files that each contain
#    one month of tick data for a single FX pair
#  - The tick directory will contain a subdirectory for each Symbol.  Each 
#    subdirectory will contain an RData file for each day containing tick data
#    for that Symbol on that day.
#  - The sec directory will look the same as the tick directory, except the
#    data will be converted to 1 second snapshots before saving.
#
# There is a lot of data here; be patient.  
#
# compressed data from May 2009 through March 2012 will use 19G of disk space
# See the very end of this script for a breakdown of disk space used.
#
## NOTE: This script requires R 2.15.0 for the paste0 function.  
## NOTE: Not intended to be used on Windows.

require(FinancialInstrument)
require(foreach)
require(doMC)

use.fasttime <- if (require(fasttime)) {
  TRUE
} else {
  message("This would be faster if the fasttime package were installed")
  message("see http://www.rforge.net/fasttime/")
  FALSE
}

################################################################################
# Begin User Parameters #
#########################
base_dir <- "~/truefx/"
## cores should not be greater than the number of gigs of available memory.
registerDoMC(cores=cores) # Can replace with a different registerDo* function
Symbols <- c("AUDJPY", "AUDNZD", "AUDUSD", "CADJPY", "CHFJPY", "EURCHF", 
             "EURGBP", "EURJPY", "EURUSD", "GBPJPY", "GBPUSD", "NZDUSD",
             "USDCAD", "USDCHF", "USDJPY")
             
## yyyymm should be a chr vector of years and months formatted YYYYMM
## The data are stored by month at truefx.com

#yyyymm <- c("201201", "201202", "201203")

curr.year <- as.numeric(format(Sys.Date(), "%Y"))
# get all the months in the current year that have already passed
curr.yyyymm <- paste0(curr.year, 
    sprintf("%02d", 1:(match(months(Sys.Date()), month.name) - 1)))
# Now combine the months from this year with the older year/months available
yyyymm <- sort(c(curr.yyyymm,
                 outer(2010:(curr.year - 1), sprintf("%02d", 1:12), paste0), 
                 paste0("2009", sprintf("%02d", 5:12))), decreasing=TRUE)
#######################
# End User Parameters #
################################################################################

if (Sys.info()[["sysname"]] == "Windows") {
  warning(paste0('This script has only been tested on linux and mac.'))
}

# if base_dir doesn't end with a forward slash, add a forward slash at the end
if (substr(base_dir, nchar(base_dir), nchar(base_dir)) != "/") {
  base_dir <- paste0(base_dir, "/")
}

# Create base_dir if it doesn't already exist as well as 3 subdirectories
dir.create(base_dir, mode="0755", showWarnings=FALSE)
dir.create(archive_dir <- paste0(base_dir, "archive/"), mode="0755", 
           showWarnings=FALSE)
dir.create(tick_dir <- paste0(base_dir, "tick/"), mode="0755", 
           showWarnings=FALSE)
dir.create(sec_dir <- paste0(base_dir, "sec/"), mode="0755", 
           showWarnings=FALSE)

# If fasttime is loaded, use fastPOSIXct, else use as.POSIXct
PosixFun <- if (use.fasttime) {
  function(x) {
    xx <- paste(paste(substring(x, 1, 4), substring(x, 5, 6), 
                      substring(x, 7, 8), sep="-"), substring(x, 10))
    fastPOSIXct(xx, "GMT")
  }
} else {
  function(x) {
    as.POSIXct(x, format="%Y%m%d %H:%M:%OS", tz="GMT")
  }
}

# set some options
if (is.null(getOption("digits.secs"))) options(digits.secs=6)
oldTZ <- Sys.getenv("TZ")
Sys.setenv(TZ='GMT') # data is stored in GMT
oldwd <- getwd()
setwd(archive_dir)

# convert "2011 09", "2011-09", or "2011/09" to "201109"
yyyymm <- substr(gsub("[ -/]", "", as.character(yyyymm)), 1, 6)

# for each of the 15 pairs, download data for each of the yyyymm months.
foreach(ym = yyyymm) %:% foreach(Symbol = Symbols) %dopar% {
  cat(Symbol, ym, "\n")
  yyyy <- substr(ym, 1, 4)
  mm <- substr(ym, 5, 6)
  zf.name <- paste0("http://www.truefx.com/dev/data/", yyyy, "/", 
                    toupper(month.name[as.numeric(mm)]), "-", yyyy, "/", 
                    Symbol, "-", yyyy, "-", mm, ".zip")
  file.create(fl <- paste0(archive_dir, Symbol, "-", yyyy, "-", mm, ".zip"))
  cat("downloading ", zf.name, "\n")
  download.file(zf.name, destfile=fl)
  cat("unzipping ", zf.name, "\n")
  uzf <- unzip(fl)
  cat("reading ", uzf, '\n')
  cat(system.time(fr <- read.csv(uzf, header=FALSE, stringsAsFactors=FALSE)), 
      "\n")
  unlink(uzf)
  id <- sub("/", "", fr[1, 1])
  cat("making index for ", id, "\n")
  idx <- PosixFun(fr[, 2])
  obj <- xts(fr[, 3:4], idx, tzone="GMT")
  #colnames(obj) <- paste(id, c("Bid.Price", "Ask.Price"), sep=".")
  colnames(obj) <- c("Bid.Price", "Ask.Price")
  tmpenv <- new.env()
  assign(id, obj, tmpenv)
  cat("saving ", id, " tick\n")
  saveSymbols.days(id, base_dir=tick_dir, extension="RData", env=tmpenv)
  assign(id, align.time(to.period(obj, "seconds", name=id, OHLC=FALSE), 1), 
         tmpenv)
  cat("saving ", id, " sec\n\n")
  saveSymbols.days(id, base_dir=sec_dir, extension="RData", env=tmpenv)
  rm("tmpenv", "obj", "fr")
}
# Restore previous settings
Sys.setenv(TZ=oldTZ)
setwd(oldwd)

## Now you should be able to
# getSymbols(Symbols, src='FI', dir=sec_dir, extension='RData', 
#            split_method='days', from='2009-05-01', days_to_omit="Saturday")
# or 
# getSymbols(Symbols, src='FI', dir=tick_dir, extension='RData', 
#            split_method='days', from='2009-05-01', days_to_omit="Saturday")



### Define instruments
#source("~/svn/blotter/pkg/FinancialInstrument/inst/parser/ISO.currencies.wiki.R")
#ccy <- unique(c(sapply(Symbols, substr, 1, 3), sapply(Symbols, substr, 4, 6)))
#
##rm_instruments()
##rm_currencies()
#define_currencies.wiki(ccy)
#fx <- exchange_rate(Symbols)
#lapply(fx, function(x) {
#  add.identifier(x, local=paste(substr(x, 1, 3), substr(x, 4, 6), sep="."))
#})
#
#lapply(fx, instrument_attr, "exchange", "TrueFX")
#lapply(fx, instrument_attr, "indexTZ", "GMT")
#lapply(fx, instrument_attr, "updated", Sys.time())
#
#saveInstruments("TrueFXinstruments.RData", dir=base_dir)
  
################################################################################
## Disk usage after downloading data for all Symbols from 200905 through 201203
#
# $ du -h truefx
# 8.7G    truefx/archive
# 147M    truefx/sec/AUDJPY
# 128M    truefx/sec/AUDNZD
# 131M    truefx/sec/AUDUSD
# 139M    truefx/sec/CADJPY
# 140M    truefx/sec/CHFJPY
# 120M    truefx/sec/EURCHF
# 143M    truefx/sec/EURGBP
# 150M    truefx/sec/EURJPY
# 156M    truefx/sec/EURUSD
# 148M    truefx/sec/GBPJPY
# 144M    truefx/sec/GBPUSD
# 77M     truefx/sec/NZDUSD
# 99M     truefx/sec/USDCAD
# 139M    truefx/sec/USDCHF
# 98M     truefx/sec/USDJPY
# 2.0G    truefx/sec
# 670M    truefx/tick/AUDJPY
# 522M    truefx/tick/AUDNZD
# 538M    truefx/tick/AUDUSD
# 626M    truefx/tick/CADJPY
# 641M    truefx/tick/CHFJPY
# 478M    truefx/tick/EURCHF
# 626M    truefx/tick/EURGBP
# 647M    truefx/tick/EURJPY
# 866M    truefx/tick/EURUSD
# 642M    truefx/tick/GBPJPY
# 574M    truefx/tick/GBPUSD
# 248M    truefx/tick/NZDUSD
# 391M    truefx/tick/USDCAD
# 609M    truefx/tick/USDCHF
# 408M    truefx/tick/USDJPY
# 8.3G    truefx/tick
# 19G     truefx

