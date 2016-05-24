###############################################
# scrape scores from Sport Affinity websites
# Portland Metro League
# url="http://oysa-2012pml.sportsaffinity.com/Tour/public/info/schedule_results2.asp?sessionguid=&flightguid={EB53B0A5-1309-4884-9A96-BE0CD7C1307C}&tournamentguid=DD593D51-4F66-4C23-9784-F00DCD2FBC91"
# Each weekends games are in a different table
# url.date.format is the format that the date is displayed on the webpage
# Friday, July 27, 2012 is date.format="%A, %B %d, %Y", but code strips off the weekname
# Scores are stored as character
###############################################
scrape.sportaffinity = function(url, file="SportAffinity", url.date.format="%B %d, %Y", date.format="%Y-%m-%d", append=FALSE, ...){
require(XML)

if(!is.character(file) | length(file)!=1)
  stop("file must be a character vector.\n",call.=FALSE)
if(!is.logical(append) | length(append)!=1)
  stop("append must be a TRUE/FALSE.\n",call.=FALSE)


tables=readHTMLTable(url, as.data.frame = TRUE, stringsAsFactors = FALSE)
tables=tables[unlist(lapply(lapply(tables,dim),length))==2]
#2012-13 tables has colnames
#game.tables = unlist(lapply(tables,function(x){ all(c("Game","Home Team","Score","Away Team") %in% colnames(x)) }))
#2013-13 they don't
#game.tables = unlist(lapply(tables,function(x){ all(c("Game","Home Team","Score","Away Team") %in% x[1,]) }))
#later in 2013-14 seem to have colnames again
game.tables = unlist(lapply(tables,function(x){ all(c("Game","Home Team","Score","Away Team") %in% x[1,]) | all(c("Game","Home Team","Score","Away Team") %in% colnames(x)) }))
games=tables[game.tables]

#Get the dates
doc = htmlParse(url)
centerNodes = getNodeSet(doc, "//center")
temp=unlist(lapply(centerNodes,xmlValue))
my.dates=sapply(temp,function(x){tmp=strsplit(x,","); str_trim(paste(tmp[[1]][2:3],collapse=", "))})
if(length(my.dates)==0){ cat("No matches on page.\n"); return() }
my.dates=as.Date(my.dates, url.date.format)
my.dates=format(my.dates, date.format)

my.table=data.frame()
for(i in 1:length(games)){
  if("Home Team" %in% games[[i]][1,]) {
    gcols=games[[i]][1,]
    grows=-1
  }else{ #probably has col names
    gcols=colnames(games[[i]]) #2012-2013
    #5-3-13; table struc has th as col names not names in row 1, don't remove first row
    grows=1:dim(games[[i]])[1]
  }
  tmp.games=games[[i]][grows,c(which(gcols=="Home Team"),which(gcols=="Home Team")+1,which(gcols=="Away Team"),which(gcols=="Away Team")+1)]
  #need to deal with PK scores which appear as 1 - 3PK
  # Problem is that entry is inconsistent
   home.score=sapply(tmp.games[,2],function(x){str_trim(str_split(x,"-")[[1]][1])})
   away.score=sapply(tmp.games[,4],function(x){str_trim(str_split(x,"-")[[1]][1])})
  #need to deal with weird characters appearing in score column
  #need to make PKs equal to NaN since score reporting is inconsistent
  home.score=sapply(home.score,function(x){paste(str_extract_all(x,"[0-9]")[[1]],collapse="")})
  away.score=sapply(away.score,function(x){paste(str_extract_all(x,"[0-9]")[[1]],collapse="")})
  tmp.table=data.frame(date=my.dates[i],home.team=tmp.games[,1],home.score=home.score,away.team=tmp.games[,3],away.score=away.score, stringsAsFactors=FALSE)
  bad.game=str_detect(tmp.games[,1],"Forfeited Game") | str_detect(tmp.games[,3],"Forfeited Game")
  tmp.table=tmp.table[!bad.game,,drop=FALSE]
  my.table=rbind(my.table,tmp.table)
}

#Set the column headings
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