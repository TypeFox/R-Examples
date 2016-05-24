################################################################
## Custom 3
## This scraper is customized to scrape the CA Coast Soccer League
##   url="http://www.coastsoccer.com/2012/LC13B.HTM"
################################################################
scrape.custom4 = function(url, file="Custom4", year=NULL, date.format="%Y-%m-%d", append=FALSE, ...){
  require(XML)

  if(missing(year)) stop("Need to know the year that the league started.\n",call.=FALSE)

  if(!is.character(file) | length(file)!=1 )
    stop("file must be a character vector.\n",call.=FALSE)
  if(!is.numeric(year) | length(year)!=1 | str_length(year)!=4)
    stop("year must be a single number and 4 digits.\n",call.=FALSE)
  if(!is.logical(append) | length(append)!=1)
    stop("append must be a TRUE/FALSE.\n",call.=FALSE)
  
tb=readHTMLTable(url,stringsAsFactors=FALSE)
tb.num=which(unlist(lapply(tb,function(x){identical(x[1,6],"FIELD")})))
  
catb=tb[[tb.num]]
catb=catb[!is.na(catb[,6]) & !str_detect(catb[,6],"FORFEITED"),c(-2,-3,-6)]
catb[,1]=sapply(catb[,1],function(x){tmp=str_remove.nonascii(x); str_proper(str_trim(str_remove(tmp,1,3))) })
catb=catb[-1,,drop=FALSE]
  months = as.numeric(format(as.Date(catb[,1],"%b %d"),"%m"))
  years=rep(year,length(months))
  if(any(diff(months)<0)) years[min(which(diff(months)<0)+1):length(months)]=year+1
  
  catb$date=paste(catb[,1],years)
  catb$date=as.Date(catb$date,"%b %d %Y")
  catb$date=format(catb$date, date.format)
catb[,2]=sapply(catb[,2],function(x){tmp=str_remove.nonascii(x); tmp=str_trim(tmp); str_replace(tmp,"\r\n","") })
catb$home.team=sapply(catb[,2],function(x){tmp=str_split(x," - ")[[1]]; tmp=str_strip.white(tmp[length(tmp)-1]); str_remove(tmp,-2)})
catb$home.score=sapply(catb[,2],function(x){tmp=str_split(x," - ")[[1]]; tmp[length(tmp)]})
catb[,3]=sapply(catb[,3],function(x){tmp=str_remove.nonascii(x); tmp=str_trim(tmp); str_replace(tmp,"\r\n","") })
catb$away.team=sapply(catb[,3],function(x){tmp=str_split(x," - ")[[1]]; tmp=str_strip.white(tmp[length(tmp)-1]); str_remove(tmp,-2)})
catb$away.score=sapply(catb[,3],function(x){tmp=str_split(x," - ")[[1]]; tmp[length(tmp)]})
catb=catb[,c(-1,-2,-3),drop=FALSE]

  extra=list(...)
  for(i in names(extra)){
    if(!(length(extra[i])==1 | length(extra[i])==dim(my.table)[1])) stop(paste("The number of values in",i,"must be 1 or equal to the number of matches."))
    my.table[i]=extra[i]
  }
  
  # Save
  if(!append) colsn=TRUE else colsn=FALSE
  if(str_sub(file, -4)!=".csv") file=paste(file,".csv",sep="")
  write.table(my.table, file=file,row.names=FALSE,col.names=colsn,append=append,sep=",",qmethod="double")
}
