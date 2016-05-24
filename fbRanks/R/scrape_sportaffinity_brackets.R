scrape.sportaffinity.brackets=function(url, file, venue=NULL, ...){
  #This is some code to scrape state cup files from affinity soccer schedule servers
  #   #EXAMPLE This is the url to the President's schedule; has the tabs with brackets
  #   url="http://2013scpdy-calsouth.affinitysoccer.com/tour/public/info/schedule_results2.asp?sessionguid=&flightguid={AA97CFBA-FB68-4E2D-BD01-DA2AF55BD819}&tournamentguid=FF92F66F-6E09-43C3-A7A7-6B052AD9C709"
  #   file="SoCal Presidents.csv"
  #   #now scrape; anything tacked on after file will be appended to scrape_affinity() arguments
  #   scrape.sportaffinity.brackets(url, file, surface="Unk",format="11v11",date.format="%m/%d/%Y")
  
  if(!str_detect(url, "schedule_results"))
    stop("Stopped in scrape.sportaffinity.brackets: url should have the words schedule_results in it.\n",.call=FALSE)
  baseurl = str_split(url,"schedule_results")[[1]][1]
  #Now scrape the page
  page <- GET(url, user_agent("httr-soccer-ranking"))
  doc = htmlParse(text_content(page))
  tableNodes = getNodeSet(doc, "//table")
  #This is the table with the bracket tabs
  the.bracket.tabs=which(unlist(lapply(tableNodes,function(x){str_sub(xmlValue(x),1,3)=="ALL"})))
  tb=tableNodes[[the.bracket.tabs]]
  #We need the html of the Node
  tb.html=as(tb,"character")
  #Now parse that html and get the a nodes
  aNodes = getNodeSet(htmlParse(tb.html), "//a")
  #The first a node is for "All" so get rid of that
  aNodes = aNodes[-1]
  #Now go through each a node and get the href and assign the bracket name
  bracket.urls=c()
  for(i in 1:length(aNodes)){
    val=xmlAttrs(aNodes[i][[1]])["href"]
    names(val)=xmlValue(aNodes[i][[1]])
    bracket.urls=c(bracket.urls,val)
  }
  
  for(i in 1:length(bracket.urls)){
    burl=paste(baseurl,bracket.urls[i],sep="")
    if(i==1){ append=FALSE }else{ append=TRUE }
    #if venue passed in, we assume the user wants venue column added
    if(!is.null(venue)){ 
      bracket=paste(venue,"Bracket",names(bracket.urls[i]))     
      scrape.sportaffinity(burl,file=file,venue=bracket, append=append, ...)
    }else{
      scrape.sportaffinity(burl,file=file, append=append, ...)
    }    
  }
}

