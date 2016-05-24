################################################################
## Custom 1
## This scraper is customized to read in json data and process
################################################################
scrape.json1 = function(url, file="Json1", date.format="%Y-%m-%d", append=FALSE,...){
  require(RJSONIO)
  require(RCurl)
  
  if(!is.character(file) | length(file)!=1 )
    stop("file must be a character vector.\n",call.=FALSE)
  if(!is.logical(append) | length(append)!=1)
    stop("append must be a TRUE/FALSE.\n",call.=FALSE)
  
  x <- getURL(url)
  x <- fromJSON(x)
  thedate = format(as.Date(x[[1]]$start), date.format)
  my.table=data.frame(date=thedate,home.team=x[[1]]$homeTeam,
                    home.score=ifelse(is.null(x[[1]]$homeScore),"NaN",x[[1]]$homeScore),
                    away.team=x[[1]]$awayTeam,
                    away.score=ifelse(is.null(x[[1]]$awayScore),"NaN",x[[1]]$awayScore))
  for(i in 2:length(x)){
    thedate = format(as.Date(x[[i]]$start), date.format)
    my.table=rbind(my.table,data.frame(date=thedate,home.team=x[[i]]$homeTeam,
        home.score=ifelse(is.null(x[[i]]$homeScore),"NaN",x[[i]]$homeScore),away.team=x[[i]]$awayTeam,
        away.score=ifelse(is.null(x[[i]]$awayScore),"NaN",x[[i]]$awayScore)))  
  }
  
  #add on extra information
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

