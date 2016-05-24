###############################################
# scrape full tables (all ages/genders) from GotSport websites
# This will dynamically determine the gender-age and put in a folder with that name
# (e.g. G01) in the basedir.
# example url="http://events.gotsport.com/events/default.aspx?EventID=23071"
# table.style 1: Most gotsoccer sites look like this one, with a game table with 8 columns
#          and header of Game, Time, Home Team, Away Team, Location
#          Date is in column 1 followed by bracket names
# table.style 2: Some gotsoccer sites have a game table with 5 columns with no header and date followed by 
#          home team names
# e.g., url="http://home.gotsoccer.com/rankings/event.aspx?EventID=24038&GroupID=236775"
# which is table.style=2, tb.num=9, url.date.format="%A, %B %d, %Y"
# Friday, July 27, 2012 is date.format="%A, %B %d, %Y"
###############################################
scrape.gotsport.main = function(url, name="test", basedir=".", tb.num=10, url.date.format="%m/%d/%Y", table.style=1, date.format="%Y-%m-%d", U12=2001, append=FALSE, ...){
  require(XML)
  require(RCurl)
  require(httr)
  
  if(!is.character(url) | length(url)!=1 )
    stop("url must be a character vector.\n",call.=FALSE)
  if(!is.numeric(tb.num) | length(tb.num)!=1)
    stop("tb.num must be a single number.\n",call.=FALSE)
  if(!is.logical(append) | length(append)!=1)
    stop("append must be a TRUE/FALSE.\n",call.=FALSE)
  
  orig.append=append
  
  #define needed subfunction
  getthelink=function(node){
    href="none" #if no link, return null
    if(!is.null(xmlChildren(node)$a)){ #if element is a link will have span element (for this the demosphere pages)
      if(!is.null(xmlChildren(node)$a)){ #that span will have an a child
        href=xmlAttrs(xmlChildren(node)$a)
        val = xmlValue(xmlChildren(node)$a)
        names(href)=val
      }
    }else{
      names(href)=str_remove.nonascii(xmlValue(xmlChildren(node)$div))
    }
    href
  }
  
  
  baseurl=str_sub(url,1,str_locate(tolower(url),"default.asp")[1]-1)
  
  #I could do just readHTMLTable(url) but I need doc later if I process the surface info
  page <- GET(url, user_agent("httr-soccer-ranking"))
  #doc = htmlParse(content(page, as = 'text'))
  doc = htmlParse(text_content(page))
  
  # read the tables and get the values in cell 1
  tb=readHTMLTable(doc, as.data.frame = TRUE, stringsAsFactors = FALSE, header=FALSE)
  lb = lapply(tb,function(x){ str_remove.nonascii(x[1,1]) })
  
  #find all the tables that have ... Groups (like Boys Groups)
  starts = c(which(str_sub(unlist(lb), -6) == "Groups"),length(lb))
  
  cat("Scraping ages.  This takes awhile....\n")
  
  for(start.i in 1:(length(starts)-1)){
    #The tables with cell 1 equal to this are the schedules
    tmp = which(lb=="FullScheduleAllStandings")
    #These are all the agetables for a group
    agetbs = tmp[tmp > starts[start.i] & tmp < starts[start.i+1]]
    
    #go through each age table
    for(tblnum in agetbs){
      #now read in the table using the function to get the urls for the schedules
      tmp=readHTMLTable(doc, elFun=getthelink, which=tblnum,as.data.frame=FALSE) #returns list
      #the urls are in the first list; it's a vector of links or "none"
      urls=tmp[[1]]
      #if the line was the name of the division, it will have value "none";
      #the name is set to the division name
      div.names=names(urls[which(urls=="none")])
      #the item right after the "none" is the link
      div.links=urls[which(urls=="none")+1]

      #Set up the file and dir name for this Gender/Age
      tmp=str_locate(div.links[1],"Gender=")[2]
      gender=unname(str_sub(div.links[1], tmp+1, tmp+1))
      tmp=str_locate(div.links[1],"Age=")[2]
      ulevel=as.numeric(unname(str_sub(div.links[1], tmp+1, tmp+2)))
      age= U12 - (ulevel-12)
      age= str_sub(as.character(age),3,4)
      genderage=paste(gender,age,sep="")
      the.dir = paste(basedir, "/", genderage,sep="")
      if(!file.exists(the.dir)) dir.create(the.dir)
      file=paste(basedir,"/",genderage,"/",name,".csv",sep="")
      
      #go through each division and scrape
      #one file for each tourny for a particular gender age
      for(l in 1:length(div.links)){
        div.link = paste(baseurl,div.links[l],sep="")
        venue=paste(name," ",div.names[l]," U",ulevel,sep="")
        
        cat(paste("Scraping",venue,genderage,"\n"))
        
        #if it is not the first division, then add onto the file
        if(l==1) append=orig.append else append=TRUE
        
        tmp=scrape.gotsport(div.link, file=file, tb.num=tb.num, url.date.format=url.date.format, table.style=table.style, date.format=date.format, append=append, venue=venue, ...)        
      }
    }   
  }
  
}