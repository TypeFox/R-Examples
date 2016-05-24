###############################################
# scrape full tables (all ages/genders) from Demosphere websites
# The age designations are variable, thus you need pass in a 3 column
# matrix with age name of div on website, U level to put in (for multi-age), div name tag to append
# e.g. BU14/15, 15, U14-15, B
# table.style=1 is whether all links are in the 'champs' table or need to get info from table below
# table.style=2 is get links from below
###############################################
scrape.demosphere.main = function(url, div.resolver, name="Demosphere", basedir=".", url.date.format="%B %d %Y", date.format="%Y-%m-%d", U12=2001, table.style=1, append=FALSE, get.surface=FALSE, ...){
  require(XML)
  require(RCurl)
  require(httr)
  
  if(!is.character(basedir) | length(basedir)!=1 )
    stop("basedir must be a character vector.\n",call.=FALSE)
  if(!is.logical(append) | length(append)!=1)
    stop("append must be a TRUE/FALSE.\n",call.=FALSE)
  if(!is.logical(get.surface) | length(get.surface)!=1)
    stop("get.surface must be a TRUE/FALSE.\n",call.=FALSE)
  if(dim(div.resolver)[2]!=4)
    stop("div.resolver must have 4 columns.\n",call.=FALSE)
  
  orig.append=append
  
  #define needed subfunctions
  getthelink=function(node){
    href="none" #if no link, return null
    if(!is.null(xmlChildren(node)$a)){ 
      if(!is.null(xmlChildren(node)$a)){ #that span will have an a child
        href=xmlAttrs(xmlChildren(node)$a)
        val = xmlValue(xmlChildren(node)$a)
        names(href)=val
      }
    }else{
      names(href)="none"
    }
    href
  }
  
  getthelink3=function(node){
    #browser()
    hrefs=c()
    href="none" #if no link, return null
    tmp=xmlChildren(node)$table
    if(!is.null(tmp)){ 
      tmp=xmlChildren(tmp)$tr
      if(!is.null(tmp)){ 
        tmp=xmlChildren(tmp)$td
        tmp1=tmp
        if(!is.null(tmp)){
          divs = which(names(xmlChildren(tmp))=="div")
          for(i in divs){
            tmp=xmlChildren(tmp1)[[i]]
            if(!is.null(tmp)){ 
              tmp=xmlChildren(tmp)$a
              if(!is.null(tmp)){ 
                href=xmlAttrs(tmp)
                val = xmlValue(tmp)
                #names(href)=val
                hrefs=c(hrefs,href,val)
              }
            }}}
      }
    }else{
      names(href)="none"
      hrefs=c(hrefs,href)
    }
    #browser()
    hrefs
  }
  
  tmp1=str_locate(url, "://")[2]
  tmp2=str_locate_all(url, "/")[[1]][,2]
  baseurl=str_sub(url,1,min(tmp2[tmp2>tmp1]))
  # read the tables and get the values in cell 1
  tb=readHTMLTable(url, as.data.frame = TRUE, stringsAsFactors = FALSE, header=FALSE)
  tblnum = which(names(tb)=="tblChamps") 
  
  #Will put all divisions for a gender-age in the same file
  used.genderages = c()
  
  cat("Scraping ages.  This takes awhile....\n")
  
  if(table.style==1){
  #now read in the table using the function to get the urls for the schedules
  tmp=readHTMLTable(url, elFun=getthelink, which=tblnum,as.data.frame=FALSE) #returns list
  #the urls are in the first list; it's a vector of links or "none"
  urls=tmp[[1]]
  urls=urls[urls!="none"]
  #if the line was the name of the division, it will have value "none";
  #the name is set to the division name
  div.names=names(urls)
  #the item right after the "none" is the link
  div.links=unname(urls)
  }else{
    #was which=tblnum+1
    tmp=readHTMLTable(url, elFun=getthelink3, which=1,as.data.frame=FALSE) #returns list
    div.names=div.links=c()
    urls=tmp[[2]]
    urls=urls[!is.na(urls) & urls!="NULL" & urls!="none"]
    for(ii in urls){
      tmp=eval(parse(text=ii))
      div.names=c(div.names,tmp[seq(2,length(tmp),2)])
      div.links=c(div.links,tmp[seq(1,length(tmp),2)])
    }
  }
  div.links=div.links[div.names!=""]
  div.names=div.names[div.names!=""]
  
    
  #Go through and scrape each age
  for(i in 1:length(div.names)){
    #Set up the file and dir name for this Gender/Age
    div.name = div.names[i]
    #Which line has the gender-age used on the website
    tmp = str_detect(div.name, div.resolver[,1])
    if(!any(tmp)) stop(paste(div.name, "is not in div.resolver.\n"))
    div.res.i = which(tmp)
    if(length(div.res.i)>1){
      tmp1=str_length(div.resolver[tmp,1])
      tmp2=which(tmp1==max(tmp1))
      div.res.i=div.res.i[tmp2]      
    }
    tmp = str_locate(div.name, div.resolver[div.res.i,1])
    div.name = str_remove(div.name, tmp[1], tmp[2])
    gender=div.resolver[div.res.i,4]
    ulevel=as.numeric(div.resolver[div.res.i,2])
    age= U12 - (ulevel-12)
    age= str_sub(as.character(age),3,4)
    genderage=paste(gender,age,sep="")
        
    the.dir = paste(basedir, "/", genderage,sep="")
    if(!file.exists(the.dir)) dir.create(the.dir)
    file=paste(basedir,"/",genderage,"/",name,".csv",sep="")
    
    div.link = paste(baseurl,div.links[i],sep="")
    #div.resolver[,3] is the age designation to use in title
    venue=str_strip.white(paste(name," ",div.name, " ", div.resolver[div.res.i,3],sep=""))
    
    cat(paste("Scraping",venue,genderage,"\n"))
    
    #if it is not the first division, then add onto the file
    tb=readHTMLTable(div.link, as.data.frame = TRUE, stringsAsFactors = FALSE, header=FALSE)
    tmp=which(names(tb)=="tblSG13c")
    tmp=as.Date(paste(unlist(tb[[tmp-1]]),"1"),"%B %Y %d")
    if(!any(is.na(tmp))){ #have months
      tmp=tmp[tmp<=Sys.Date()]
      tmp=format(tmp,"%Y%m")
      div.link=paste(str_remove(div.link,-4),tmp,".html",sep="")
    }
    for(div.linki in div.link){ #deal with multiple div.links due to months
      if(genderage %in% used.genderages) append=TRUE else append=orig.append
      tmp=try(scrape.demosphere(div.linki, file=file, url.date.format=url.date.format, date.format=date.format, append=append, get.surface=get.surface, venue=venue, ...)  )      
    if(class(tmp)=="try-error"){
      cat("There was an error scraping ", div.linki, "\n")
      cat("Call was: scrape.demosphere(div.linki, file=file, url.date.format=url.date.format, date.format=date.format, append=append, get.surface=get.surface, venue=venue, ...)\n")
      browser()
    }
      used.genderages = c(used.genderages, genderage) 
    }
    #update the list of used gender ages
  }
}   
