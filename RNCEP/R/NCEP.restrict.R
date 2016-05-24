NCEP.restrict <-
function(wx.data, years2remove=NULL, months2remove=NULL, days2remove=NULL, hours2remove=NULL, other2remove=NULL, set2na=TRUE){

## wx.data three dimensional array as returned by gather.wx.data()

## hours2remove can be any combination of c(0,6,12,18)
## months2remove should be given using numerical naming convention i.e. c(10,12)
## years2remove should be given in YYYY format i.e. c(2001,2002)
## other2remove is for specific combinations of year,month,day, and hour given -- must be given in the format ("%Y-%m-%d %H") !!!!

## set2na is a logical. Should the values be set to NA (default) or removed completely?
## If set2na is something other than TRUE or FALSE, replace observations with value of set2na.

if(set2na == TRUE) { 

for(i in 1:length(unlist(dimnames(wx.data)[3]))){
	if(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][1] != "XXXX"){
		if(any(years2remove == as.numeric(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][1]))){
		wx.data[,,i] <- NA }}
	if(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][2] != "XX"){
		if(any(months2remove == as.numeric(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][2]))){
		wx.data[,,i] <- NA }}
	if(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][3] != "XX"){
		if(any(days2remove == as.numeric(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][3]))){
		wx.data[,,i] <- NA }}
	if(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][4] != "XX"){
		if(any(hours2remove == as.numeric(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][4]))){
		wx.data[,,i] <- NA }}
	if(any(gsub(x=other2remove, pattern = "-", replacement = " ") == gsub(x=dimnames(wx.data)[[3]][i], pattern = "_", replacement = " "))){
		wx.data[,,i] <- NA }
	}
} else


if(set2na == FALSE) {
removal <- c()
	
for(i in 1:length(unlist(dimnames(wx.data)[3]))){
	if(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][1] != "XXXX"){
		if(any(years2remove == as.numeric(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][1]))){
		removal <- append(removal, i) }}
	if(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][2] != "XX"){
		if(any(months2remove == as.numeric(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][2]))){
		removal <- append(removal, i) }}
	if(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][3] != "XX"){
		if(any(days2remove == as.numeric(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][3]))){
		removal <- append(removal, i) }}
	if(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][4] != "XX"){
		if(any(hours2remove == as.numeric(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][4]))){
		removal <- append(removal, i) }}
	if(any(gsub(x=other2remove, pattern = "-", replacement = " ") == gsub(x=dimnames(wx.data)[[3]][i], pattern = "_", replacement = " "))){
		removal <- append(removal, i) }
}
	wx.data <- wx.data[,,-removal]
} else
	

for(i in 1:length(unlist(dimnames(wx.data)[3]))){
	if(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][1] != "XXXX"){
		if(any(years2remove == as.numeric(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][1]))){
		wx.data[,,i] <- set2na }}
	if(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][2] != "XX"){
		if(any(months2remove == as.numeric(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][2]))){
		wx.data[,,i] <- set2na }}
	if(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][3] != "XX"){
		if(any(days2remove == as.numeric(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][3]))){
		wx.data[,,i] <- set2na }}
	if(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][4] != "XX"){
		if(any(hours2remove == as.numeric(strsplit(unlist(dimnames(wx.data)[3])[i], split="_")[[1]][4]))){
		wx.data[,,i] <- set2na }}
	if(any(gsub(x=other2remove, pattern = "-", replacement = " ") == gsub(x=dimnames(wx.data)[[3]][i], pattern = "_", replacement = " "))){
		wx.data[,,i] <- set2na }
	}

return(wx.data)
}

