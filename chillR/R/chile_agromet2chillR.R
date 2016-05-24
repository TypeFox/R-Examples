#' Convert a weather file downloaded from the Chilean Agromet website to chillR
#' format
#' 
#' Convert downloaded weather data into a data frame that makes running other
#' chillR functions easy.
#' 
#' Processing the data with this function will make the data work well with the
#' remainder of this package.
#' 
#' @param downloaded_weather_file full path of a weather file downloaded from
#' the Chilean Agromet website (http://agromet.inia.cl/) as an alleged Excel
#' file (it has some formatting issues).
#' @param drop_most boolean variable indicating if most columns should be
#' dropped from the file. If set to TRUE (default), only essential columns for
#' running chillR functions are retained.
#' @return a data.frame with weather data, according to the downloaded file
#' provided as input. If drop_most is FALSE, all columns from the original
#' dataset are preserved, although some column names are adjusted to chillR's
#' preferences ("Year","Month","Day","Tmin","Tmax","Tmean","Prec", if these
#' columns are present). If drop_most is TRUE, only columns likely to be of
#' interest to chillR users are retained.
#' @note Many databases have data quality flags, which may sometimes indicate
#' that data aren't reliable. These are not considered by this function!
#' @author Eike Luedeling
#' @references The chillR package:
#' 
#' Luedeling E, Kunz A and Blanke M, 2013. Identification of chilling and heat
#' requirements of cherry trees - a statistical approach. International Journal
#' of Biometeorology 57,679-689.
#' @keywords utilities
#' @examples
#' 
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2005),])  # this line is
#' #only here to make the example run, even without downloading a file
#' 
#' # FOLLOW THE INSTRUCTIONS IN THE LINE BELOW THIS; AND THEN  RUN THE LINE AFTER THAT (without the #)
#' # download an Excel file from the website and save it to disk (path: {X}) 
#' #weather<-fix_weather(chile_agromet2chillR({x}))
#' 
#' hourtemps<-stack_hourly_temps(weather, latitude=50.4)
#' chilling(hourtemps,305,60)
#' 
#' @export chile_agromet2chillR
chile_agromet2chillR<-function(downloaded_weather_file,drop_most=TRUE)
{xml <- htmlParse(downloaded_weather_file)
 ns <- getNodeSet(xml, '//tr')
 ns<-ns[-c(1:5)]
 days<-lapply(ns, function(x) {event <- t(xmlToDataFrame(downloaded_weather_file))})
 caption<-days[[1]]
 caption<-caption[which(!is.na(caption))]
 days<-days[which(sapply(days,length)==as.numeric(names(sort(table(sapply(days,length)),decreasing=TRUE)[1])))]
 res<-data.frame(matrix(unlist(days),ncol=length(days[[2]]),byrow=T))
 comma_replace<-function(y) {st<-strsplit(as.character(y),",")[[1]]
 if(length(st)==2) y<-as.numeric(st[1])+as.numeric(st[2])/10
 if(st[[1]]=="-") y<-NA
 return(y)}
 for(i in 1:ncol(res)) res[,i]<-unlist(sapply(res[,i],comma_replace))
 for(i in 1:ncol(res)) if(!length(strsplit(as.character(res[1,i]),"-")[[1]])>1) res[,i]<-as.numeric(res[,i])
 colnames(res)<-caption
 res[,"Year"]<-as.numeric(sapply(strsplit(as.character(res$FECHA),"-"),function(x) x[3]))
 res[,"Month"]<-as.numeric(sapply(strsplit(as.character(res$FECHA),"-"),function(x) x[2]))
 res[,"Day"]<-as.numeric(sapply(strsplit(as.character(res$FECHA),"-"),function(x) x[1]))
 colnames(res)[grep("xim",colnames(res))]<-"Tmax"
 colnames(res)[grep("nim",colnames(res))]<-"Tmin"
 colnames(res)[grep("Aire",colnames(res))]<-"Tmean"
 colnames(res)[grep("Prec",colnames(res))]<-"Prec"
 colnames(res)[grep("FECHA",colnames(res))]<-"Date"
 name_seq<-c("Year","Month","Day","Tmin","Tmean","Tmax","Prec")
 namess<-name_seq[which(name_seq %in% colnames(res))]
 othernames<-colnames(res)[which(!colnames(res) %in% namess)]
 if(drop_most) res<-res[,namess] else res<-res[,c(namess,othernames)]
return(res)}


