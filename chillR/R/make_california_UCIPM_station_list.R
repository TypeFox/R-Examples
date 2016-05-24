#' Makes a list of the UC IPM weather stations
#' 
#' Makes a list of the weather stations contained in the UC IPM database, with
#' geographic coordinates. This requires parsing through quite a few websites,
#' because the coordinates don't seem to be stored in one central (and easily
#' accessible) place. Hence this is much slower than one might expect. A
#' shortcut is the california_stations dataset supplied with chillR, which
#' contains the result of running this function in February 2016. The default
#' in the other relevant functions will be the use of this pre-stored list, but
#' if the current station coverage is needed, this function can help. Having
#' said this, station coverage probably won't change very rapidly, so in most
#' cases, the california_stations dataset should be enough.
#' 
#' 
#' @return a data.frame containing stations from the California UC IPM database
#' (), with the following columns: "Name", "Code", "Interval", "Lat", "Long",
#' "Elev".
#' @author Eike Luedeling
#' @references The chillR package:
#' 
#' Luedeling E, Kunz A and Blanke M, 2013. Identification of chilling and heat
#' requirements of cherry trees - a statistical approach. International Journal
#' of Biometeorology 57,679-689.
#' @keywords utilities
#' @examples
#' 
#' #cali_stats<-make_california_UCIPM_station_list()
#' @export make_california_UCIPM_station_list
make_california_UCIPM_station_list<-function()  
  {docu<-htmlParse("http://ipm.ucdavis.edu/WEATHER/wxactstnames.html")
  els = getNodeSet(docu, "//body//table")[[2]]
  els = getNodeSet(els, "//table//tr")
  nores=TRUE
  for(i in 1:length(els))
  {x<-xmlToDataFrame(els[i])
  if(length(x)==3)
  {colnames(x)<-c("Name","Code","Interval")
  if(nores) {res<-x;nores=FALSE} else res<-rbind(res,x)}}
  
  for(l in 1:nrow(res))
  {docu<-htmlParse(paste("http://ipm.ucdavis.edu/calludt.cgi/WXSTATIONDATA?STN=",res$Code[l],sep=""))
  els = getNodeSet(docu, "//table")[[2]]
  els = getNodeSet(els, "//tr")[[6]]
  positionstring<-getChildrenStrings(els)[1]
  suppressWarnings(sp<-as.numeric(strsplit(positionstring," ")$td))
  sp<-sp[which(!is.na(sp))]
  res[l,"Lat"]<-sp[1]+sp[2]/60
  res[l,"Long"]<-sp[3]+sp[4]/60
  if(length(grep("min W",positionstring))>0) res[l,"Long"]<-(-res[l,"Long"])
  if(length(grep("min S",positionstring))>0) res[l,"Lat"]<-(-res[l,"Lat"]) 
  res[l,"Elev"]<-as.numeric(strsplit(as.character(getChildrenStrings(els)[3])," ")[[1]][2])*0.3048}
  return(res)
}

