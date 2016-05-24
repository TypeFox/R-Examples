get_cfsr_latlon<-function(declat,declon,emailaddr,timeoff=0,interppow=2){
#
# Grabs historical CFSR data through time for a given lat and lon (over unfrozen land surface) using the service at: drfuka.org
#
  options(timeout=600)
  url=paste("http://www.cfsr.tamu-cornell.drfuka.org/swat-cfsr-v03.pl?lat=",declat,"&lon=",declon,"&timeoff=",timeoff,"&interppow=",interppow,"&.submit=Submit",sep="")
  urlline=grep("data/data",readLines(url),value=T)
  urlgz=strsplit(urlline,"\"")[[1]][2]

  download.file(urlgz,"junk.gz")
  hist_wx=read.csv(zz <- gzfile("junk.gz"),header=T,colClasses=c("character","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"))
  hist_wx$date=as.Date(hist_wx$date,format="%Y-%m-%d")
  file.remove("junk.gz")

  colnames(hist_wx)=c("DATE","TMX","TMN","PRECIP","WIND","AVGRH","MAXRH","MINRH","SOLAR")

  return(hist_wx)
}

