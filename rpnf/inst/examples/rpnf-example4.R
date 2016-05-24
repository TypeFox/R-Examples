#' This example explains how to create Point & Figure Charts for all commodities of a given composite index, e.g. the GDAXI
#' Moreover relative strength and bullish percent charts are created, too.
#' 
### Initialize libraries
library(rpnf) # Load rpnf library
library(quantmod) # Load quantmod library
stockData <- new.env() #Make a new environment for quantmod to store data in

### Define wrapper for quantmod download function
downloadDatas <- function(symbols=c("GOOG")) {
  startDate = as.Date(Sys.Date()-10*365) #Specify period of time we are interested in
  endDate = as.Date(Sys.Date()+1) # End date is today
  getSymbols(symbols, env = stockData, src = "yahoo", from = startDate, to = endDate)   #Download the stock history (for all tickers)
  ### rename xts-columns
  for (s in symbols) 
    eval(parse(text=paste("names(stockData$",sub("^\\^","",s),")<-c(\"open\",\"high\",\"low\",\"close\",\"volume\",\"adjusted\")",sep="")))
}

### Define wrapper for getting an appropriate symbol table
getSymbolTable <- function(index) {
  yahoo.url <- "http://finance.yahoo.com/d/quotes.csv?s="
  download.file(url=paste(yahoo.url,"@",index,"&f=snbaopl1&e=.csv",sep=""),destfile="temp.csv")
  # read temp.csv file
  symbol <- read.csv("temp.csv", header=F)[,1]
  # return result
  data.frame(index,symbol)
}

# Function to show progress on 
printProgress <- function(i,n,step=1,start.time=NULL) {
  if (i %% step == 0 | i == n) {
    percentDone <- i/n;
    if (!is.null(start.time)) {
      timeElapsed <- as.numeric(difftime(Sys.time(), start.time, units="secs"));
      totalTimePrediction <- timeElapsed/percentDone;
      timeToFinish <- totalTimePrediction-timeElapsed;
      linesPerSecond <- i/timeElapsed;
      cat(paste(i, "/", n," \t(", round(100*percentDone, digits=2), "%) \t",round(linesPerSecond,digits=2)," lines/s\t",round(totalTimePrediction,digits=2),"s ETT \t",round(timeToFinish,digits=2),"s ETF\n", sep=""));
    } else {
      cat(paste(i, "/", n," \t (", round(100*percentDone, digits=2), "%)"))       
    }
  }
}

### The example code starts
boxsize <- getLogBoxsize(percent=2.5)
log <- TRUE
# Define (yahoo) index symbol (with "^") to be processed
index <- "^GDAXI" # e.g. GDAXI, DJI, see http://finance.yahoo.com for more
# download stock quotes for index
downloadDatas(c(index))
index.xts <- eval(parse(text=paste("OHLC(stockData$",sub("^\\^","",index),")",sep="")))
# create PNF information for index
indexPnf <- pnfprocessor(index.xts[,2],index.xts[,3],index(index.xts),boxsize=boxsize,log=log)
# Get list of composite symbols
symbolTable <- getSymbolTable(index)
# download stock quotes for composite symbols
downloadDatas(as.character(symbolTable$symbol))
# generate point and figure information for composite index
symbolPnf <- data.frame()
symbolRS <- data.frame()
rs <- data.frame()
i<-0
start.time <- Sys.time()
for (symbol in symbolTable[,2]) {
  print(paste("processing symbol: ",symbol))
  #data <- symbolQuotes[symbolQuotes$symbol==symbol,]
  #symbolPnf <- rbind(symbolPnf, cbind(symbol,pnfprocessor(data$high,data$low,date=as.Date(data$date),boxsize=boxsize,log=log)))
  symbol.xts <- eval(parse(text=paste("OHLC(stockData$",sub("^\\^","",symbol),")",sep="")))
  symbolPnf <- rbind(symbolPnf, cbind(symbol,
                                      pnfprocessor(symbol.xts[,2],symbol.xts[,3],date=index(symbol.xts),boxsize=boxsize,log=log)))
  # determine relative strength
  rs.xts <- symbol.xts/index.xts
  ### use new calcRS function
  symbolRS <- rbind(symbolRS, cbind(symbol,
                                    pnfprocessor(rs.xts[,2],rs.xts[,3],date=index(rs.xts),boxsize=boxsize,log=log,style="rs")))
  i<-i+1
  printProgress(i,n=length(symbolTable[,2]),step=1,start.time)
}

# generate bullish percent chart
boxsizeBP <- 2
bptable<-table(symbolPnf$date,symbolPnf$status.bs)
bp<-100*bptable[,"Buy"]/(bptable[,"Buy"]+bptable[,"Sell"])
bpPnf <- pnfprocessor(high=bp,date=as.Date(names(bp)),boxsize=boxsizeBP,log=FALSE,style="bp")
# generate ascending percent chart
asctable<-table(symbolPnf$date,symbolPnf$status.xo)
asc<-100*asctable[,"X"]/(asctable[,"X"]+asctable[,"O"])
ascPnf <- pnfprocessor(high=asc,date=as.Date(names(asc)),boxsize=boxsizeBP,log=FALSE,style="bp")
# generate plots
i<-0
start.time <- Sys.time()
for (s in symbolTable[,2]) {
  print(paste("plotting symbol: ",s))
  #oldpar <- par()
  # create plotting canvas
  png(filename=paste(s,".png"),width=3200,height=1800)
  par(mfrow=c(4,1))
  # plot symbol chart
  pnfplot(symbolPnf[symbolPnf$symbol==s,],boxsize=boxsize,log=log,main=paste("Commodity chart: ",s))
  # plot symbol rs
  pnfplot(symbolRS[symbolRS$symbol==s,],boxsize=boxsize,log=log,main=paste("Relative strength chart: ",s," vs. ",index))
  # plot index bp
  pnfplot(bpPnf,boxsize=boxsizeBP,log=F,main=paste("Bullish percent chart: ",index))
  # plot index asc
  pnfplot(ascPnf,boxsize=boxsizeBP,log=F,main=paste("Ascending percent chart: ",index))
  # restore plot settings
  dev.off()
  #par(oldpar)
  i<-i+1
  printProgress(i,n=length(symbolTable[,2]),step=1,start.time)
}

