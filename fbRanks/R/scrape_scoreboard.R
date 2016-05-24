#Select the html of the table from a ScoreBoard page
#Save to a file
scrape.scoreboard=function(html.file, file="ScoreBoard", url.date.format="%a %m/%d/%Y", date.format="%Y-%m-%d", append=FALSE, get.surface=FALSE, ...){

  if(!is.character(file) | length(file)!=1 )
    stop("file must be a character vector.\n",call.=FALSE)
  if(!is.logical(append) | length(append)!=1)
    stop("append must be a TRUE/FALSE.\n",call.=FALSE)  
  
  test=readHTMLTable(html.file)
  my.table=c()
  surface=c()
  for(tbl in names(test)[-1]){
    tmp=test[[tbl]]
    colnames(tmp)=apply(tmp[1,],2,as.character)
    tmp=tmp[-1,]
    surf=rep("Unk",dim(tmp)[1])
    surf[str_detect(tmp$Venue,"Turf")]="Turf"
    surf[str_detect(tmp$Venue,"Grass")]="Grass"
    tmp=tmp[,c("Date","Home Team","Score","AwayTeam"),drop=FALSE]
    tmp[,5]=unlist(lapply(str_split(tmp$Score,"-"),function(x){x[2]}))
    tmp[,3]=unlist(lapply(str_split(tmp$Score,"-"),function(x){x[1]}))
    surface=c(surface,surf)
    tmp[,1]=format(as.Date(tmp$Date, url.date.format), date.format)
    my.table=rbind(my.table,tmp)
  }
  colnames(my.table)=c("date","home.team","home.score", "away.team", "away.score")
  
  extra=list(...)
  for(i in names(extra)){
    if(!(length(extra[i])==1 | length(extra[i])==dim(my.table)[1])) stop(paste("The number of values in",i,"must be 1 or equal to the number of matches."))
    my.table[i]=extra[i]
  }
  if(get.surface)
    my.table$surface=surface
  
  # Save
  if(!append) colsn=TRUE else colsn=FALSE
  if(str_sub(file, -4)!=".csv") file=paste(file,".csv",sep="")
  write.table(my.table, file=file,row.names=FALSE,col.names=colsn,append=append,sep=",",qmethod="double")
}