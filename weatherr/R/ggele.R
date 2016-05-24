ggele = function(lat=0,lon=0, output=c('elevation','elevation/resolution','all')) {
output <- match.arg(output)
if ((length(lat) != length(lon)) | (any(abs(lat)>90)) | (any(abs(lon)>180))) stop('Longitude and latitude should have equal length and within the valid range')
inpar = paste(paste0(lat,',', lon),collapse = ' | ')
url = paste0("http://maps.googleapis.com/maps/api/elevation/json?locations=", inpar)
u = paste0(paste0(readLines(url), collapse = "\n"), "\n")
tmp = try(RJSONIO::fromJSON(u))
if (class(tmp) == 'try-error') stop('Error: fail to obtain data from Google Elevation API')
if (tmp$status != 'OK') stop('Request denied')
switch(output, 'elevation' = {y=sapply(tmp$results,function(x) x$elevation);names(y)=1:length(tmp$results);y}, 'elevation/resolution' = {y=sapply(tmp$results,function(x) c(x$elevation,x$resolution));y=as.data.frame(y);colnames(y)=1:length(tmp$results);rownames(y) = c('elevation','resolution');y}, all = tmp)
}