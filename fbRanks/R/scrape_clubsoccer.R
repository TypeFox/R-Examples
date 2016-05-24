###############################################
# scrape scores from US Club Soccer 
# url="http://clubsoccer.us/TTSchedules.aspx?tid=USCLUBG&year=2013&stid=USCLUBG&syear=2012&div=U12M01"
###############################################
scrape.usclub = function(url, file="USClub", url.date.format="%A%m/%d/%Y", date.format="%Y-%m-%d", append=FALSE, ...){
require(XML)

if(!is.character(file) | length(file)!=1 )
  stop("file must be a character vector.\n",call.=FALSE)
if(!is.logical(append) | length(append)!=1)
  stop("append must be a TRUE/FALSE.\n",call.=FALSE)

tb=readHTMLTable(url, as.data.frame = TRUE, stringsAsFactors = FALSE)
tb.start=which(names(tb)=="ctl00_ContentPlaceHolder1_TableBtn")+3
#tb.end=which(names(tb)=="ctl00_ContentPlaceHolder1_ChangesTable")-1
tb.end=tb.start+min(which(!(names(tb)[tb.start:length(tb)]=="NULL")))-2
game.tbs=tb[seq(tb.start,tb.end,2)]
for(i in 1:length(game.tbs)){
  this.tb=game.tbs[[i]]
  this.tb=this.tb[!is.na(this.tb[,"V2"]),,drop=FALSE]
  the.dates=as.Date(this.tb[seq(1,dim(this.tb)[1],2),"V2"], url.date.format)
  the.dates=format(the.dates, date.format)
  the.home.team=this.tb[seq(1,dim(this.tb)[1],2),"V7"]
  the.home.score=this.tb[seq(1,dim(this.tb)[1],2),"V5"]
  the.away.team=this.tb[seq(2,dim(this.tb)[1],2),"V3"]
  the.away.score=this.tb[seq(2,dim(this.tb)[1],2),"V1"]
  this.table=data.frame(date=the.dates,home.team=the.home.team, home.score=the.home.score, away.team=the.away.team, away.score=the.away.score)
  if(i==1) my.table=this.table
  else my.table=rbind(my.table,this.table)
}

# Set the column headings
colnames(my.table)=c("date","home.team","home.score", "away.team", "away.score")


#Replace missing scores with NaN
my.table$home.score[my.table$home.score==""]=NaN
my.table$away.score[my.table$away.score==""]=NaN

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