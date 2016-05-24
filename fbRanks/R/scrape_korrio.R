################################################################
## Scrape scores off sites run off Korrio
################################################################
scrape.korrio=function(url, file="Korrio", url.date.format="%B %Y %a %d", date.format="%Y-%m-%d", append=FALSE, get.surface=FALSE, ...){
require(httr)
require(XML)

if(!is.character(file) | length(file)!=1 )
  stop("file must be a character vector.\n",call.=FALSE)
if(!is.logical(append) | length(append)!=1)
  stop("append must be a TRUE/FALSE.\n",call.=FALSE)

page <- GET(url, user_agent("httr-soccer-ranking"))
my.table <- readHTMLTable(text_content(page), as.data.frame=TRUE, stringsAsFactors=FALSE)
if(length(my.table)==0){
  cat(text_content(page))
  stop("scrape.korrio did not return a valid table list")
}
my.table = my.table[[2]]
my.table = my.table[,c(1,4,5,6),drop=FALSE]

is.date=(1:dim(my.table)[1])[is.na(my.table$Home)]
for(i in 1:(length(is.date)-1)){
if(length((is.date[length(is.date)]+1):(is.date[i+1]-1))!=0)
  my.table[(is.date[i]+1):(is.date[i+1]-1),1]=paste(my.table[is.date[i],1],my.table[(is.date[i]+1):(is.date[i+1]-1),1])
}
if(length(is.date[length(is.date)]:dim(my.table)[1])!=0)
  my.table[is.date[length(is.date)]:dim(my.table)[1],1]=paste(my.table[is.date[length(is.date)],1],my.table[is.date[length(is.date)]:dim(my.table)[1],1])
#get rid of the rows with month
my.table=my.table[-is.date,,drop=FALSE]
my.table[,1]=as.Date(my.table[,1], url.date.format) #read in the dates
my.table[,1]=format(my.table[,1], date.format) #reformat the dates

# Fix the score column
my.table[,3]=as.character(my.table[,3])
my.table[,5]=NaN
scr=lapply(strsplit(my.table[,3],"-"),str_trim)
for(i in 1:length(scr)){
  if(identical(scr[[i]],"TBD") | length(scr[[i]])!=2){
  my.table[i,3]=NaN
}else{
  my.table[i,5]=scr[[i]][2]
  my.table[i,3]=scr[[i]][1]
}
if(any(str_detect(c(my.table[i,3],my.table[i,5]),"[*]"))){ #not official
    my.table[i,c(3,5)]=NaN
  }
}

#Fix the column headings
colnames(my.table)=c("date","home.team","home.score", "away.team", "away.score")

#set the extra columns of information
extra=list(...)
for(i in names(extra)){
  if(!(length(extra[i])==1 | length(extra[i])==dim(my.table)[1])) stop(paste("The number of values in",i,"must be 1 or equal to the number of matches."))
  my.table[i]=extra[i]
}


# Set surface if requested; takes time
if(get.surface){
  cat("Getting surface information.  This takes awhile....\n")
  page <- GET(url, user_agent("httr-soccer-ranking"))
  doc = htmlParse(text_content(page))
  aNodes = getNodeSet(doc, "//a")
  links.to.details=which(unlist(lapply(aNodes,xmlValue))=="G")
  details.urls=c()
  for(i in links.to.details){
    details.urls=c(details.urls,xmlAttrs(aNodes[i][[1]])["href"])
  }
  surface=c()
  for(url in details.urls){
    page <- GET(url, user_agent("httr-soccer-ranking"))
    doc = htmlParse(content(page, as = 'text'))
    fieldTypeNodes = getNodeSet(doc, "//div[@class='span2 pull-right cal_event_detail_fieldtype']")    
    #tmp=xmlValue(getNodeSet(x,"//strong")[[8]]) #works but below also works
    tmp=str_split(str_strip.white(xmlValue(fieldTypeNodes[[1]]))," ")[[1]][2]
    surf="Unk"
    if(tmp=="Tu") surf="Turf"
    if(tmp=="Gr") surf="Grass"
    surface=c(surface,surf)
  }
  my.table$surface=surface
}

# Save
if(!append) colsn=TRUE else colsn=FALSE
if(str_sub(file, -4)!=".csv") file=paste(file,".csv",sep="")
write.table(my.table, file=file,row.names=FALSE,col.names=colsn,append=append,sep=",",qmethod="double")
}