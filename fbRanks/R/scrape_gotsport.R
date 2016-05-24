###############################################
# scrape scores from GotSport websites
# url="http://events.gotsport.com/events/schedule.aspx?EventID=23071&Gender=Boys&Age=12"
# table.style 1: Most gotsoccer sites look like this one, with a game table with 8 columns
#          and header of Game, Time, Home Team, Away Team, Location
#          Date is in column 1 followed by bracket names
# table.style 2: Some gotsoccer sites have a game table with 5 columns with no header and date followed by 
#          home team names
# e.g., url="http://home.gotsoccer.com/rankings/event.aspx?EventID=24038&GroupID=236775"
# which is table.style=2, tb.num=9, url.date.format="%A, %B %d, %Y"
# Friday, July 27, 2012 is date.format="%A, %B %d, %Y"
###############################################
scrape.gotsport = function(url, file="GotSport", tb.num=10, url.date.format="%m/%d/%Y", table.style=1, date.format="%Y-%m-%d", append=FALSE, ...){
require(XML)

if(!is.character(file) | length(file)!=1 )
  stop("file must be a character vector.\n",call.=FALSE)
if(!is.numeric(tb.num) | length(tb.num)!=1)
  stop("tb.num must be a single number.\n",call.=FALSE)
if(!is.logical(append) | length(append)!=1)
  stop("append must be a TRUE/FALSE.\n",call.=FALSE)


if(table.style==1){ tb=readHTMLTable(url, as.data.frame = TRUE, stringsAsFactors = FALSE, header=FALSE)
}else{ tb=readHTMLTable(url, as.data.frame = TRUE, stringsAsFactors = FALSE, header=FALSE) }
if(length(tb)<tb.num){
  cat(paste(file,": The selected tb.num=",tb.num,"is bigger than the number of tables.\nThe list of tables is being sent to the output.\nIf the game table has 5 columns and no header, try table.style=2.\n"))
  return(tb)
}
if(table.style==1){
if(dim(tb[[tb.num]])[2]!=8){
  ok.tables=which(unlist(lapply(tb,function(x){ifelse(is.null(dim(x)), FALSE, dim(x)[2]==8)})))
  cat(paste(file,": The selected tb.num=",tb.num,"does not look right.  It should have 8 columns.\n",
             ifelse(length(ok.tables)==0,"No tables have 8 columns. Try table.style=2.\n",
                 paste("tb.num",ok.tables,"have 8 columns.\n",collapse=" ")),
                    "The list of tables is being sent to the output.\nIf the game table has 5 columns and no header, try table.style=2.\n"))
  return(tb)
}
my.table=tb[[tb.num]][,c(-2,-7,-8),drop=FALSE]
}
if(table.style==2){
  if(dim(tb[[tb.num]])[2]!=5){
    cat(paste("The selected tb.num=",tb.num,"does not look right.  It should have 5 columns and no header.\nPass in a different table number (try 12 or 10).\nThe list of tables is being sent to the output.\nIf the game table has 8 columns and header, try table.style=1."))
    return(tb)
  }
  my.table=tb[[tb.num]][,c(-3),drop=FALSE]
  my.table=cbind(my.table[,1],my.table) #add date column since it is missing
}
# Fix the date column
my.table[,1]=as.Date(my.table[,1], url.date.format) #read in the date
my.table[,1]=format(my.table[,1], date.format) #reformat the date
is.date=(1:dim(my.table)[1])[!is.na(my.table[,1])]
if(length(is.date)==0) stop("There's a problem.  No dates in table.\n")
#This for loop is to put the dates in column 1
for(i in 1:length(is.date)){ #for all date rows
  last.date.row=is.date[length(is.date)]
  if((i+1)>length(is.date)){ #we are at the last date
    next.date.row=dim(my.table)[1]+1 #then fill in to the last row of table
  }else{ next.date.row=is.date[i+1] }
  this.date.row=is.date[i]
  #fill in from row after this.date.row to next.date.row-1 (or end of table if at last date)
  my.table[(this.date.row+1):(next.date.row-1),1]=my.table[this.date.row,1]
}

# Set the column headings
colnames(my.table)=c("date","home.team","home.score", "away.team", "away.score")

#Clean up NAs
my.table=my.table[!is.na(my.table$away.score),,drop=FALSE]
my.table=my.table[!is.na(my.table$date),,drop=FALSE]
my.table=my.table[!(my.table$home.team=="Home Team"),,drop=FALSE]

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