################################################################
## Scrape scores off sites by Demosphere
## table style 1 http://www.pugetsoundpremierleague.com/schedules/Spring2013/57712924.html
## table style 2 http://www.norcalpremier.com/schedules/59258181/59258272.html
################################################################
scrape.demosphere = function(url, file="Demosphere", url.date.format="%B %d %Y", date.format="%Y-%m-%d", table.style=1, year=NULL, append=FALSE, get.surface=FALSE, ...){
  require(XML)
  require(RCurl)
  require(httr)
  
  if(!is.character(file) | length(file)!=1 )
    stop("file must be a character vector.\n",call.=FALSE)
  if(!is.logical(append) | length(append)!=1)
    stop("append must be a TRUE/FALSE.\n",call.=FALSE)
  if(!is.logical(get.surface) | length(get.surface)!=1)
    stop("get.surface must be a TRUE/FALSE.\n",call.=FALSE)
  if(!(table.style %in% c(1,2)))
    stop("table.style must be a 1 or 2.\n",call.=FALSE)
  if(table.style==2 & is.null(year))
    stop("if table.style is 2, you must specify the year as a 4 digit number.\n",call.=FALSE)
  if(table.style==2 & !is.null(year)){
    if(!is.numeric(year)) stop("if table.style is 2, you must specify the year as a 4 digit number.\n",call.=FALSE)
    if(str_length(as.character(year))!=4) stop("if table.style is 2, you must specify the year as a 4 digit number.\n",call.=FALSE)
  }
  #I could do just readHTMLTable(url) but I need doc later if I process the surface info
  #   if(get.surface){
  #     page <- GET(url, user_agent("httr-soccer-ranking"))
  #   #doc = htmlParse(content(page, as = 'text'))
  #   doc = htmlParse(text_content(page))
  #   }else{ doc=url }
  
  doc=url
  
  # read the tables and select one labeled tbListGames2
  tables=readHTMLTable(doc)
  # the match info we want is in columns 1,4,5,6
  if(table.style==1){
    my.table=tables$tblListGames2[,c(1,4,5,6),drop=FALSE]
    # get rid of rows that are not matches.  Those will have GAME# in column 1
    my.table=my.table[my.table[,1]!="GAME#",]
    #set up the away.score column
    my.table[,5]=as.character(NaN)  
    #Set the column headings
    colnames(my.table)=c("date","home.team","home.score", "away.team", "away.score")    
    # Set up the date column
    my.table$date=as.character(my.table$date)
    #date column sometimes has funky characters in the front
    my.table$date=sapply(my.table$date,function(x){tmp=strsplit(x,", ");paste(unlist(tmp)[-1],collapse=" ")})
    
    
    my.table$date=as.Date(my.table$date, url.date.format) #read in the dates
    my.table$date=format(my.table$date, date.format) #reformat the dates
    #now figure out which rows are dates and which are not; non-dates will be NAs
    is.date=(1:dim(my.table)[1])[!is.na(my.table$date)]
    if(length(is.date)>1){ #more than one date
      for(i in 1:(length(is.date)-1)){
        if(length((is.date[length(is.date)]+1):(is.date[i+1]-1))!=0)
          my.table[(is.date[i]+1):(is.date[i+1]-1),"date"]=my.table[is.date[i],"date"]
      }
    }
    if(length(is.date[length(is.date)]:dim(my.table)[1])!=0) #there are games for last date
      my.table[is.date[length(is.date)]:dim(my.table)[1],"date"]=my.table[is.date[length(is.date)],"date"]
    my.table=my.table[-is.date,,drop=FALSE]
    
    # Get rid of any rows that are NAs; those were associated with the match date and not a match result
    my.table=my.table[!is.na(my.table$away.team),,drop=FALSE]
  }#end of table.style 1
  if(table.style==2){
    my.table=tables$tblListGames2[,c(4,6,7,8),drop=FALSE]
    # get rid of rows that are not matches.  Those will have NA in column 1
    my.table=my.table[!is.na(my.table[,1]),]
    #set up the away.score column
    my.table[,5]=as.character(NaN)  
    #Set the column headings
    colnames(my.table)=c("date","home.team","home.score", "away.team", "away.score")
    #Set up the dates
    my.table$date=as.character(my.table$date)
    my.table$date=str_remove.nonascii(my.table$date)
    my.table=my.table[!is.na(as.Date(my.table$date,"%b %d")),]
    my.table$date=paste(my.table$date,year)
    my.table$date=as.Date(my.table$date, "%b %d %Y") #read in the dates
    my.table$date=format(my.table$date, date.format) #reformat the dates 
    # Get rid of any date rows that are NAs; those are not a match result
    my.table=my.table[!is.na(my.table$date),,drop=FALSE]    
  }
  
  # Fix the score column
  my.table$home.score=as.character(my.table$home.score)
  scr=strsplit(my.table$home.score,"-")
  for(i in 1:length(scr)){
    #FT or FFT means forfeit
    if(identical(scr[[i]],"vs") | str_detect(paste(scr[[i]],collapse=" "),"FT")){
      my.table$home.score[i]=NaN
    }else{
      my.table$away.score[i]=scr[[i]][2]
      my.table$home.score[i]=scr[[i]][1]
    }  
  }
  
  extra=list(...)
  for(i in names(extra)){
    if(!(length(extra[i])==1 | length(extra[i])==dim(my.table)[1])) stop(paste("The number of values in",i,"must be 1 or equal to the number of matches."))
    my.table[i]=extra[i]
  }
  
  if(get.surface){
    cat("Getting surface information.  This takes awhile....\n")
    #above we read in all the tables at the match url and put it in tables
    #the table we want is tblListGames2; figure out which table that is
    tblnum=match("tblListGames2",names(tables))
    #read in the match table again but this time as a list so as.data.frame doesn't monkey with it
    testvals=readHTMLTable(doc, which=tblnum,as.data.frame=FALSE) #returns list
    #get the locations which are in column 7 but since didn't convert to dataframe, it's in element 7
    locationnames=testvals[[7]]
    #The values in column 7 when row is not a game will be Location or NA; get rid of those
    locationnames=locationnames[locationnames!="Location" & !is.na(locationnames)]
    #now get the link associated with each location.  The location is hotlinked and we need the url
    #to do this we need to write a function for readHTMLTable to use for each elmement
    getthelink=function(node){
      val=NULL #if no link, return null
      tmp=node
      if(!is.null(xmlChildren(tmp)$div)) tmp=xmlChildren(tmp)$div
      if(!is.null(xmlChildren(tmp)$div)) tmp=xmlChildren(tmp)$div
      if(!is.null(xmlChildren(tmp)$span)) tmp=xmlChildren(tmp)$span
      tmp=xmlChildren(tmp)$a
      if(!is.null(tmp)) #that span will have an a child
        val=xmlAttrs(tmp) #href is at attribute of the "a" child
      val
    }
    #now read in the table using the function to get the urls for the fields
    test=readHTMLTable(doc, elFun=getthelink, which=tblnum,as.data.frame=FALSE) #returns list
    #the urls are not in element 7
    details.urls=test[[7]]
    #if the line was not a game, will be NA or "NULL";
    details.urls=details.urls[!is.na(details.urls) & !(details.urls=="NULL")]
    
    #now that we have the field urls and field names, we can load each of the field pages and look up the surface information
    surface=c()
    for(i in 1:length(details.urls)){
      url=details.urls[i]
      fieldname=str_strip.white(locationnames[i])
      #if the fieldname is "" then it means it is not listed and surface is "Unk"
      if(fieldname==""){
        surface=c(surface,"Unk")
      }else{ #the fieldname is listed so we can try to get the surface information
        #the http is script.  This will run the script and get the resulting html
        x <- getURL(url)
        #it returns text.  We need to parse this into html with a node tree
        doc = htmlParse(x)
        #now we can read in the tables.
        fieldtable=readHTMLTable(doc) #this will be a dataframe
        #The table we want iunfortunately it doesn't have a name
        fieldtable.num=max(which(unlist(lapply(fieldtable,function(x){"FIELD NAME" %in% x[,1]}))))
        fieldtable=fieldtable[[fieldtable.num]]
        #the fieldname and surface are in columns 1 and 4 respectively; first 2 lines are header info
        first.row = which(fieldtable[,1]=="FIELD NAME")
        fieldtable=fieldtable[first.row:dim(fieldtable)[1],,drop=FALSE]
        if(!any(fieldtable[1,]=="SURFACE")) fieldtable=cbind(fieldtable,X=c("SURFACE",rep("Unk",dim(fieldtable)[1]-1)))
        needed.cols=c(which(fieldtable[1,]=="FIELD NAME"),which(fieldtable[1,]=="SURFACE"))
        #remove the header and non-needed columns
        fieldtable=fieldtable[-1,needed.cols]
        #find the line where the field name matches that in the match page; strip off any extra white space
        
        fieldline=match(fieldname,str_strip.white(fieldtable[,1]))
        #if there happens not to be a match, don't crash
        if(!is.na(fieldline)){
          #get the surface information
          tmp=str_strip.white(as.character(fieldtable[fieldline,2]))
          if(tmp=="-") tmp="Unk"
          tmp=str_proper(tmp) #str_proper is just converting to having the first letter capitalized
          #Guard against stuff like "Field Turf" instead of just Turf
          if(str_detect(tmp,"Turf") | str_detect(tmp,"Truf")) tmp="Turf"
          if(str_detect(tmp,"Grass")) tmp="Grass"
          if(!(tmp %in% c("Turf","Grass","Unk"))) tmp="Unk"
          surface=c(surface,str_proper(tmp)) 
        }else{ surface=c(surface,"Unk") }
      }
    }
    my.table$surface=surface
  }
  
  #Finally check that no team names==""
  my.table=my.table[my.table$home.team!="",,drop=FALSE]
  
  # Save
  if(!append) colsn=TRUE else colsn=FALSE
  if(str_sub(file, -4)!=".csv") file=paste(file,".csv",sep="")
  write.table(my.table, file=file,row.names=FALSE,col.names=colsn,append=append,sep=",",qmethod="double",na="NaN")
}