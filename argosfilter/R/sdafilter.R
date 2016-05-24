.packageName<-"argosfilter"

`sdafilter` <-
function(lat, lon, dtime, lc, vmax=2,ang=c(15,25),distlim=c(2500,5000))
{
row_id<-1:length(lat);row_id
original_data<-data.frame(row_id,lat,lon,dtime,lc);original_data[1:10,]
#length(original_data$lat)
#------------------------------------------------------------------------
# calculate vmask
#------------------------------------------------------------------------
mfilter<-vmask(lat,lon,dtime,vmax)
mfilter
#length(which(mfilter=="removed"))
original_data<-cbind(original_data,mfilter)
#original_data[1:10,]
#length(original_data$lat)

#------------------------------------------------------------------------
# remove Z positions and vmask="removed" (if dist>5000m to the prev pos)
#------------------------------------------------------------------------
dist<-dist.prev(lat,lon)
original_data<-cbind(original_data,dist);#original_data
rem<-which((mfilter=="removed" & dist>5000) | lc==-9 | lc=="z" | lc=="Z" )
if(length(rem)!=0) filt_data<-original_data[-rem,] else filt_data<-original_data
#filt_data[1:10,]
#length(original_data$lat)
#length(filt_data$lat)

#------------------------------------------------------------------------
# remove angles<=angles in the list & dist>distances in the list
#------------------------------------------------------------------------

remx<-c(1,2)
while (!is.na(remx[1])){
	latx<-filt_data$lat
	lonx<-filt_data$lon
	#length(latx)

	# angle to the previous and next position
	angs<-internal.angles(latx,lonx)
	#angs[1:10]

	# distance to the next position
	dnext<-dist.next(latx,lonx)
	#dnext[1:10]

	# distance to the previous position
	dprev<-dist.prev(latx,lonx)
	#dprev[1:10]

	rem_string=character(1);rem_string
	for (i in 1:length(ang) ){
		string<-paste("(angs<=",ang[i]," & dprev>",distlim[i]," & dnext>",distlim[i],")",sep="")
		if (i==1) rem_string<-string else rem_string=paste(rem_string,"|",string)
	}
	rem_string
	eval(parse(text=rem_string))
	remx<-which( eval(parse(text=rem_string)) )
	#remx
	if (!is.na(remx[1])) filt_data<-filt_data[-remx,]
	#filt_data[15:35,]
	#length(filt_data$lat)
}
#filt_data[1:10,]
#length(original_data$lat)
#length(filt_data$lat)

#------------------------------------------------------------------------
# Find removed rows
#------------------------------------------------------------------------
not_rem<-match(filt_data$row_id,original_data$row_id);#not_rem
removed<-original_data[not_rem*-1,]
#filt_data[1:20,]
#removed[1:20,]
sda_filter<-character(length(original_data$lat));#sda_filter
sda_filter[not_rem*-1]<-"removed";#sda_filter
sda_filter[not_rem]<-"not";#sda_filter
extremes<-c(1,2,length(original_data$lat)-1,length(original_data$lat));#extremes
sda_filter[extremes]<-"end_location";#sda_filter
sda_filter

}

