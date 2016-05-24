###############################################
# scrape full tables (all ages/genders) from SportAffinity websites with no brackets
#Example
#url=http://wys-2013rcl.affinitysoccer.com/tour/public/info/accepted_list.asp?tournamentguid={4AAB9406-8A3D-4069-9CC3-C85037E1ABC6}
###############################################
# url.date.format is the format that the date is displayed on the webpage
# Friday, July 27, 2012 is date.format="%A, %B %d, %Y", but code strips off the weekname
# Scores are stored as character
# U.designation is the text right before the age;  so Boys Under 11  U.designation = "Under " or "der "; BU11 is "U"
# name.delimiter and name.skip help form the division name
# Say divisions are like Boys Under 11 Foobar and you want to use Foobar as the division; name.delimiter="Under ", name.skip=3
# if you want 11 Foobar, name.skip=0
# if you want Under 11 Foobar, name.delimiter="Boys" name.skip=1
###############################################
scrape.sportaffinity.main = function(url, name="SportAffinity", basedir=".", url.date.format="%B %d, %Y", date.format="%Y-%m-%d", U12=2001, append=FALSE, add.to.base.url="tour/public/info/", ..., U.designation="Under ", name.delimiter="Under ", name.skip=3){
  require(XML)
  require(RCurl)
  require(httr)
  
  if(!is.character(basedir) | length(basedir)!=1 )
    stop("basedir must be a character vector.\n",call.=FALSE)
  if(!is.logical(append) | length(append)!=1)
    stop("append must be a TRUE/FALSE.\n",call.=FALSE)
  
  orig.append=append
  
  getthelink3=function(node){
     hrefs=c()
    href="none" #if no link, return null
    tmp=xmlChildren(node)$div
    if(!is.null(tmp)){
      tmp=xmlChildren(tmp)$table
      if(!is.null(tmp)){ 
        trs = which(names(xmlChildren(tmp))=="tr")
        tmp1=tmp
        for(i in trs){
            tmp=xmlChildren(tmp1)[[i]]
              tds = which(names(xmlChildren(tmp))=="td")
              if(length(tds)!=0){
                if(!is.null(xmlChildren(tmp[[tds[3]]])$b)){
                if(!is.null(xmlChildren(xmlChildren(tmp[[tds[3]]])$b)$a)){
                val=xmlValue(tmp[[tds[1]]])
                href=xmlAttrs(xmlChildren(xmlChildren(tmp[[tds[3]]])$b)$a)
                hrefs=c(hrefs,href,val)
                }
              }}
            }
        }
      }
    hrefs
  }
  
  tmp1=str_locate(url, "://")[2]
  tmp2=str_locate_all(url, "/")[[1]][,2]
  baseurl=str_sub(url,1,min(tmp2[tmp2>tmp1]))
  baseurl=paste(baseurl,add.to.base.url,sep="")
  
  
  #Will put all divisions for a gender-age in the same file
  used.genderages = c()
  
  cat("Scraping ages.  This takes awhile....\n")
  
    tmp=readHTMLTable(url, elFun=getthelink3, which=8, as.data.frame=FALSE) #returns list
    for(i in 1:length(tmp)){
      tmp[[i]]=tmp[[i]][!is.na(tmp[[i]]) & !(tmp[[i]]=="NULL")]
    }
    tmp=unlist(tmp)
    div.names=tmp[seq(2,length(tmp),2)]
    div.links=tmp[seq(1,length(tmp),2)] 
    
  #Go through and scrape each age
  for(i in 1:length(div.names)){
    #Set up the file and dir name for this Gender/Age
    div.name = div.names[i]
    gender=str_sub(div.name,1,1)
    if(str_detect(div.name, "HS") | str_detect(div.name, "High School")){
      ulevel=18
    }else{
    ulevel=str_sub(div.name,str_locate(div.name,U.designation)[2]+1,str_locate(div.name,U.designation)[2]+2)
    ulevel=as.numeric(ulevel)
    }
    age= U12 - (ulevel-12)
    age= str_sub(as.character(age),3,4)
    genderage=paste(gender,age,sep="")
    
    div.name = str_remove(div.name, 1, str_locate(div.name,name.delimiter)[2]+name.skip)
    
        
    the.dir = paste(basedir, "/", genderage,sep="")
    if(!file.exists(the.dir)) dir.create(the.dir)
    file=paste(basedir,"/",genderage,"/",name,".csv",sep="")
    
    div.link = paste(baseurl,div.links[i],sep="")
    #div.resolver[,3] is the age designation to use in title
    venue=str_strip.white(paste(name," ",genderage," ",div.name,sep=""))
    
    cat(paste("Scraping",venue,"\n"))
    
    #if it is not the first division, then add onto the file
    if(genderage %in% used.genderages) append=TRUE else append=orig.append
    
    tmp=try(scrape.sportaffinity(div.link, file=file, url.date.format=url.date.format, date.format=date.format, append=append, venue=venue, ...))
    if(class(tmp)=="try-error"){
      cat("There was an error scraping ", div.link, "\n")
      browser()
    }
    #update the list of used gender ages
    used.genderages = c(used.genderages, genderage) 
  }
}   
