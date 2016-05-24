locationforecast <- function(lat,lon,elevation=NULL,location=NULL,exact=TRUE,tz=Sys.timezone()) {
if (!is.null(location)) {
latlon = as.numeric(rev(geocode(location=location,source = "google")))
if (any(is.null(latlon))) stop('Error: no match location')
lat = latlon[1]; lon = latlon[2]
elevation = fromJSON(paste0('http://maps.googleapis.com/maps/api/elevation/json?locations=',lat,',',lon,'&sensor=false'))$results[[1]]$elevation
}
if (any(!is.numeric(c(lat,lon,elevation)))) stop('Error: lat, lon and elevation have to be numeric')
lat = lat[1]; lon = lon[1]; elevation=elevation[1]
msl = ifelse(is.null(elevation),'',paste0(';msl=',round(elevation,digits=0)))
url = paste0('http://api.met.no/weatherapi/locationforecast/1.9/?lat=',lat,';lon=',lon,msl)
x =  paste0(paste0(readLines(url), collapse = "\n"), "\n")
x = xmlRoot(xmlParse(x))
timefrom = with_tz(strptime(xpathSApply(x,'//time/@from'),'%Y-%m-%dT%H:%M:%S',tz='GMT'),tzone=tz)
timeto = with_tz(strptime(xpathSApply(x,'//time/@to'),'%Y-%m-%dT%H:%M:%S',tz='GMT'),tzone=tz)
tdiff = difftime(timeto,timefrom,units='hours')
if (exact) {
v = which(tdiff==0)
alltime = xpathSApply(x,'//time')
allexact = sapply(v,function(y) alltime[[y]])
temp = sapply(allexact,function(z) list(point=xmlToList(z[[1]])))
temperature = as.numeric(sapply(1:length(temp),function(x) {x1 = temp[[x]]$temperature['value']; ifelse(!is.null(x1),x1,NA)}))
windDirection = as.numeric(sapply(1:length(temp),function(x) {x1 = temp[[x]]$windDirection['deg']; ifelse(!is.null(x1),x1,NA)}))
windSpeed_mps = as.numeric(sapply(1:length(temp),function(x) {x1 = temp[[x]]$windSpeed['mps']; ifelse(!is.null(x1),x1,NA)}))
windSpeed_beaufort = as.numeric(sapply(1:length(temp),function(x) {x1 = temp[[x]]$windSpeed['beaufort']; ifelse(!is.null(x1),x1,NA)}))
windSpeed_name = sapply(1:length(temp),function(x) {x1 = temp[[x]]$windSpeed['name']; ifelse(!is.null(x1),x1,NA)})
windGust = as.numeric(sapply(1:length(temp),function(x) {x1 = temp[[x]]$windGust['mps']; ifelse(!is.null(x1),x1,NA)}))
humidity = as.numeric(sapply(1:length(temp),function(x) {x1 = temp[[x]]$humidity['value']; ifelse(!is.null(x1),x1,NA)}))
pressure = as.numeric(sapply(1:length(temp),function(x) {x1 = temp[[x]]$pressure['value']; ifelse(!is.null(x1),x1,NA)}))
cloudiness = as.numeric(sapply(1:length(temp),function(x) {x1 = temp[[x]]$cloudiness['percent']; ifelse(!is.null(x1),x1,NA)}))
lowClouds = as.numeric(sapply(1:length(temp),function(x) {x1 = temp[[x]]$lowClouds['percent']; ifelse(!is.null(x1),x1,NA)}))
mediumClouds = as.numeric(sapply(1:length(temp),function(x) {x1 = temp[[x]]$mediumClouds['percent']; ifelse(!is.null(x1),x1,NA)}))
highClouds = as.numeric(sapply(1:length(temp),function(x) {x1 = temp[[x]]$highClouds['percent']; ifelse(!is.null(x1),x1,NA)}))
dewpointTemperature = as.numeric(sapply(1:length(temp),function(x) {x1 = temp[[x]]$dewpointTemperature['value']; ifelse(!is.null(x1),x1,NA)}))
data.frame(time=timefrom[v], temperature=temperature,windDirection=windDirection, windSpeed_mps=windSpeed_mps, windSpeed_beaufort=windSpeed_beaufort,windSpeed_name=windSpeed_name, windGust=windGust, humidity=humidity, pressure=pressure, cloudiness=cloudiness, lowClouds=lowClouds, mediumClouds=mediumClouds, highClouds=highClouds, dewpointTemperature=dewpointTemperature)
} else
{
v = which(tdiff!=0)
alltime = xpathSApply(x,'//time')
allexact = sapply(v,function(y) alltime[[y]])
temp = sapply(allexact,function(z) list(point=xmlToList(z[[1]])))
precipitation = as.numeric(sapply(1:length(temp),function(x) {x1 = temp[[x]]$precipitation['value']; ifelse(!is.null(x1),x1,NA)}))
minTemperature = as.numeric(sapply(1:length(temp),function(x) {x1 = temp[[x]]$minTemperature['value']; ifelse(!is.null(x1),x1,NA)}))
maxTemperature = as.numeric(sapply(1:length(temp),function(x) {x1 = temp[[x]]$maxTemperature['value']; ifelse(!is.null(x1),x1,NA)}))
weather_id = sapply(1:length(temp),function(x) {x1 = temp[[x]]$symbol['id']; ifelse(!is.null(x1),x1,NA)})
data.frame(timefrom=timefrom[v], timeto=timeto[v], interval = tdiff[v], precipitation=precipitation, minTemperature=minTemperature, maxTemperature=maxTemperature, weather_id=weather_id, stringsAsFactors=F)
}
}
