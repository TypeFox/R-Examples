################################################################
## Custom 2
## This scraper is customized to scrape the Vancouver Metro League
## http://www.bcsoccercentral.com/Schedules/MSL/dirg6.htm
################################################################
scrape.custom2 = function(url, file="Custom2", year=NULL, date.format="%Y-%m-%d", append=FALSE, ...){
  require(XML)
  
  if(missing(year)) stop("Need to know the year that the league started.\n",call.=FALSE)
  
  if(!is.character(file) | length(file)!=1 )
    stop("file must be a character vector.\n",call.=FALSE)
  if(!is.numeric(year) | length(year)!=1 | str_length(year)!=4)
      stop("year must be a single number and 4 digits.\n",call.=FALSE)
  if(!is.logical(append) | length(append)!=1)
    stop("append must be a TRUE/FALSE.\n",call.=FALSE)
  
  tb = readHTMLTable(url,stringsAsFactors=FALSE,colClasses="character")[[1]]
  tb = tb[,-2]
  asc <- iconv(tb$Team, "latin1", "ASCII")
  team.rows <- !(is.na(asc) | asc != tb$Team)
  date.char = str_proper(str_strip.white(names(tb)))
  months = as.numeric(format(as.Date(date.char,"%b %d"),"%m"))
  years=rep(year,length(months))
  if(any(diff(months)<0)) years[min(which(diff(months)<0)+1):length(months)]=year+1
  date.char=paste(date.char,years)
  dates = as.Date(date.char,"%b %d %Y")
  date.cols=which(!is.na(dates))
  team.names=tb$Team[team.rows]
  my.table=c()
  for(row in which(team.rows)){
    tb[row,]=str_strip.white(as.character(tb[row,]))
    scores=str_remove.nonascii(as.character(tb[row+1,]))
    home.team=tb$Team[[row]]
    for(col in date.cols){
      if(str_sub(tb[row,col],1,1)=="H"){ #home game
        away.team=tb$Team[away.team=tb[,1]==str_sub(tb[row,col],2)]        
        date=dates[[col]]
        home.score=away.score=NaN
        if(str_detect(scores[[col]],"-")){
          tmp=str_split(scores[[col]],"-")[[1]]
          home.score=str_strip.white(tmp[1])
          away.score=str_strip.white(tmp[2])
        }
        my.table=rbind(my.table,data.frame(date=date,home.team=home.team,
                                           home.score=home.score,away.team=away.team,away.score=away.score,stringsAsFactors=FALSE))  
      }
    }
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

