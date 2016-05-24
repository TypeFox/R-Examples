################################################################
## Custom 3
## This scraper is customized to scrape the USYSA Regionals
## http://tournaments.usyouthsoccer.org/events/2012-Region-IV-Championships/Schedule/U12-Boys/Division+1/
################################################################
scrape.custom3 = function(url, file="Custom3", year=NULL, date.format="%Y-%m-%d", append=FALSE, ...){
  require(XML)
  
  if(missing(year)) stop("Need to know the year that the league started.\n",call.=FALSE)

  if(!is.character(file) | length(file)!=1 )
    stop("file must be a character vector.\n",call.=FALSE)
  if(!is.numeric(year) | length(year)!=1 | str_length(year)!=4)
    stop("year must be a single number and 4 digits.\n",call.=FALSE)
  if(!is.logical(append) | length(append)!=1)
    stop("append must be a TRUE/FALSE.\n",call.=FALSE)
  
  tb = readHTMLTable(url,stringsAsFactors=FALSE,colClasses="character")
  match.tbls=unlist(lapply(tb,function(x){any(x[1,]=="Score")}))
  
  my.table=c()
  for(i in which(match.tbls)){
    tbl=tb[[i]]
    colnames(tbl)=as.character(tbl[1,])
    tbl=tbl[-1,,drop=FALSE]
    date=paste(tbl$Date,year)
    date=as.Date(date,"%a %b %d %Y")
    home.team=tbl$Home
    away.team=tbl$Away
    scores=str_strip.white(str_remove.nonascii(tbl$Score,sub=" "))
    home.score=sapply(scores,function(x){str_split(x," ")[[1]][1]})
    away.score=sapply(scores,function(x){str_split(x," ")[[1]][2]})
    my.table=rbind(my.table,data.frame(date=date,home.team=home.team,
                                       home.score=home.score,away.team=away.team,away.score=away.score,stringsAsFactors=FALSE))      
  }
  my.table$date=format(my.table$date, date.format) #reformat the dates
      
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

