################################################################
## Custom 1
## This scraper is customized to deal with score table html that is malformed
## and is missing <tr> tags
## first.td.tag is the # of the first <td> tag associated with the table.
## last.td.tag is the last <td> tag in the table relative to all <td> tags in page; 0 means last tag; 2 means 2 before last tag
## Also score don't have dates on this site, just numbers for weeks; weeks should be in the date.format desired.
## Found on http://www.district-3.org
## Found on North Puget Sound League; 2012-2013 first.td.tag=4; 2013-2014 first.td.tag=3
################################################################
scrape.custom1 = function(url, file="Custom1", weeks=NULL, first.td.tag=3,last.td.tag=7,td.per.row=5, append=FALSE,...){
require(XML)

if(!is.character(file) | length(file)!=1)
  stop("file must be a character vector.\n",call.=FALSE)
for(elem in c("first.td.tag","last.td.tag","td.per.row")){
  el = get(elem)
  if(!is.numeric(el) | length(el)!=1)
  stop(paste(elem,"must be a single number.\n"),call.=FALSE)
}
if(!is.logical(append) | length(append)!=1)
  stop("append must be a TRUE/FALSE.\n",call.=FALSE)

if(missing(weeks)) stop("This type of webpage has week numbers instead of dates.\n  You need to specify what date each week number corresponds to.\n")
dates=weeks

doc = htmlParse(url)

#Normally we would just do this:
#tableNodes=getNodeSet(doc, "//table")
#tb = readHTMLTable(tableNodes[[2]])

#The following hack gets around this
tableNodes = getNodeSet(doc, "//td")

tb.rows = seq(first.td.tag,(length(tableNodes)-last.td.tag),5)
db.rows = 1:length(tb.rows)
my.table = data.frame(date=rep(1,length(db.rows)),home.team="a",
home.score=1,away.team="a",away.score=1,stringsAsFactors=FALSE)

#set date
for(i in db.rows){
val=as.numeric(xmlValue(tableNodes[[tb.rows[i]]]))
my.table$date[i]=dates[val]
}

#set home team name
tb.rows = seq((first.td.tag+2),(length(tableNodes)-last.td.tag),5)
for(i in db.rows){
val=as.character(xmlValue(tableNodes[[tb.rows[i]]]))
val=str_sub(val,end=-3)
my.table$home.team[i]=val
}

#set away team name
tb.rows = seq((first.td.tag+4),(length(tableNodes)-last.td.tag),5)
for(i in db.rows){
val=as.character(xmlValue(tableNodes[[tb.rows[i]]]))
val=str_sub(val,end=-1)
my.table$away.team[i]=val
}

#set home score
tb.rows = seq((first.td.tag+1),(length(tableNodes)-last.td.tag),5)
for(i in db.rows){
val=as.numeric(xmlValue(tableNodes[[tb.rows[i]]]))
my.table$home.score[i]=val
}

#set away score
tb.rows = seq((first.td.tag+3),(length(tableNodes)-last.td.tag),5)
for(i in db.rows){
val=as.numeric(xmlValue(tableNodes[[tb.rows[i]]]))
my.table$away.score[i]=val
}

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

