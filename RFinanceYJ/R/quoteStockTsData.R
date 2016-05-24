#API
quoteStockTsData <- function(x, since=NULL,start.num=0,date.end=NULL,time.interval='daily')
{
  time.interval <- substr(time.interval,1,1)
  function.stock <- function(quote.table.item){
    if( xmlSize(quote.table.item) < 5) return(NULL) 
    d <- convertToDate(xmlValue(quote.table.item[[1]]),time.interval)
    o <- as.number(xmlValue(quote.table.item[[2]]))
    h <- as.number(xmlValue(quote.table.item[[3]]))
    l <- as.number(xmlValue(quote.table.item[[4]]))
    c <- as.number(xmlValue(quote.table.item[[5]]))
    v <- ifelse(xmlSize(quote.table.item) >= 6,as.number(xmlValue(quote.table.item[[6]])),0)
    a <- ifelse(xmlSize(quote.table.item) >= 7,as.number(xmlValue(quote.table.item[[7]])),0)
    return(data.frame(date=d,open=o,high=h,low=l,close=c,volume=v, adj_close=a))
  }
  return(quoteTsData(x,function.stock,since,start.num,date.end,time.interval,type="stock"))
}
quoteFundTsData <- function(x, since=NULL,start.num=0,date.end=NULL,time.interval='daily')
{
  time.interval <- substr(time.interval,1,1)
  function.fund <- function(quote.table.item){
    d <- convertToDate(xmlValue(quote.table.item[[1]]),time.interval)
    if(time.interval=='monthly'){
      d <- endOfMonth(d)
    }
    c <- as.number(xmlValue(quote.table.item[[2]]))
    v <- as.number(xmlValue(quote.table.item[[3]]))
    return(data.frame(date=d,constant.value=c,NAV=v))
  }
  return(quoteTsData(x,function.fund,since,start.num,date.end,time.interval,type="fund"))
}
quoteFXTsData <- function(x, since=NULL,start.num=0,date.end=NULL,time.interval='daily')
{
  time.interval <- substr(time.interval,1,1)
  function.fx <- function(quote.table.item){
    d <- convertToDate(xmlValue(quote.table.item[[1]]),time.interval)
    o <- as.number(xmlValue(quote.table.item[[2]]))
    h <- as.number(xmlValue(quote.table.item[[3]]))
    l <- as.number(xmlValue(quote.table.item[[4]]))
    c <- as.number(xmlValue(quote.table.item[[5]]))
    return(data.frame(date=d,open=o,high=h,low=l,close=c))
  }
  return(quoteTsData(x,function.fx,since,start.num,date.end,time.interval,type="fx"))
}
######  private functions  #####
#get time series data from Yahoo! Finance.
quoteTsData <- function(x,function.financialproduct,since,start.num,date.end,time.interval,type="stock"){
  r <- NULL
  result.num <- 51
  financial.data <- data.frame(NULL)
  #start <- (gsub("([0-9]{4,4})-([0-9]{2,2})-([0-9]{2,2})","&c=\\1&a=\\2&b=\\3",since))
  #end   <- (gsub("([0-9]{4,4})-([0-9]{2,2})-([0-9]{2,2})","&f=\\1&d=\\2&e=\\3",date.end))
  start <- (gsub("([0-9]{4,4})-([0-9]{2,2})-([0-9]{2,2})","&sy=\\1&sm=\\2&sd=\\3",since))
  end   <- (gsub("([0-9]{4,4})-([0-9]{2,2})-([0-9]{2,2})","&ey=\\1&em=\\2&ed=\\3",date.end))

  if(!any(time.interval==c('d','w','m'))) stop("Invalid time.interval value")
  
  extractQuoteTable <- function(r,type){
    if(type %in% c("fund","fx")){
      tbl <- r[[2]][[2]][[7]][[3]][[3]][[9]][[2]]
    }
    else{
      tbl <- r[[2]][[2]][[7]][[3]][[3]][[10]][[2]]
    }
    return(tbl)
  }
  
  while( result.num >= 51 ){
    start.num <- start.num + 1
    quote.table <- NULL
    quote.url <- paste('http://info.finance.yahoo.co.jp/history/?code=',x,start,end,'&p=',start.num,'&tm=',substr(time.interval,1,1),sep="")
  
    try( r <- xmlRoot(htmlTreeParse(quote.url,error=xmlErrorCumulator(immediate=F))), TRUE)
    if( is.null(r) ) stop(paste("Can not access :", quote.url))

    #try( quote.table <- r[[2]][[1]][[1]][[16]][[1]][[1]][[1]][[4]][[1]][[1]][[1]], TRUE )
    try( quote.table <- extractQuoteTable(r,type), TRUE )
    
    if( is.null(quote.table) ){
      if( is.null(financial.data) ){
        stop(paste("Can not quote :", x))
      }else{
         financial.data <- financial.data[order(financial.data$date),]
         return(financial.data)
      }
    }

    size <- xmlSize(quote.table)
    for(i in 2:size){
      financial.data <- rbind(financial.data,function.financialproduct(quote.table[[i]]))
    }
    
    result.num <- xmlSize(quote.table)
    Sys.sleep(1)
  }
  financial.data <- financial.data[order(financial.data$date),]
  return(financial.data)  
}
#convert string formart date to POSIXct object
convertToDate <- function(date.string,time.interval)
{
  #data format is different between monthly and dialy or weekly
  if(any(time.interval==c('d','w'))){
    result <- gsub("^([0-9]{4})([^0-9]+)([0-9]{1,2})([^0-9]+)([0-9]{1,2})([^0-9]+)","\\1-\\3-\\5",date.string)
  }else if(time.interval=='m'){
    result <- gsub("^([0-9]{4})([^0-9]+)([0-9]{1,2})([^0-9]+)","\\1-\\3-01",date.string)
  }
  return(as.POSIXct(result))
}
#convert string to number.
as.number <- function(string)
{
  return(as.double(as.character(gsub("[^0-9.]", "",string))))
}
#return end of month date.
endOfMonth <- function(date.obj)
{
  startOfMonth     <- as.Date(format(date.obj,"%Y%m01"),"%Y%m%d")
  startOfNextMonth <- as.Date(format(startOfMonth+31,"%Y%m01"),"%Y%m%d")
  return(startOfNextMonth-1)
}



