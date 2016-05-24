.packageName<-"argosfilter"

vmask<-function(lat, lon, dtime, vmax)
{
row_id<-1:length(lat)
v<-numeric(length(lat))
dset<-data.frame(row_id,lat,lon,dtime,v,row.names = NULL)

dset2<-dset
n_int=0
maxi=100

ascending=TRUE
curr_peak=0
curr_null=0
n_peaks=0
pos_peak=0


while(maxi>vmax){
	n_int=n_int+1 # number of iteractions
	lat<-dset2$lat
	lon<-dset2$lon
	dtime<-dset2$dtime

	#-----------------------------------------------------
	#calculate velocities v[i]
	for (i in 3:(length(lat)-2)) {
		v_2=distance_m(lat[i],lon[i],lat[i-2],lon[i-2])/as.numeric(difftime(dtime[i],dtime[i-2],units = "secs")+1)
		v_1=distance_m(lat[i],lon[i],lat[i-1],lon[i-1])/as.numeric(difftime(dtime[i],dtime[i-1],units = "secs")+1)
		v1=distance_m(lat[i],lon[i],lat[i+1],lon[i+1])/as.numeric(difftime(dtime[i+1],dtime[i],units = "secs")+1)
		v2=distance_m(lat[i],lon[i],lat[i+2],lon[i+2])/as.numeric(difftime(dtime[i+2],dtime[i],units = "secs")+1)
		v[i]=sqrt(sum(v_2^2, v_1^2, v1^2, v2^2)/4)
		dset2$v[i]=v[i]
	}

	#-----------------------------------------------------
	# get the peaks in v[i]
	ascending=TRUE
	curr_peak=0
	curr_null=0
	n_peaks=0
	pos_peak=0
	for (i in 3:(length(lat)-2)) {
		if (ascending) {                                
			if (v[i]>curr_peak) curr_peak=v[i]	else {
				ascending = FALSE;              
				curr_null = v[i];
				pos_peak=cbind(pos_peak,i-1);
				n_peaks=n_peaks+1
			}

		} else {
			if (v[i] < curr_null) curr_null = v[i] else {
				ascending = TRUE # previous point was a minimum, now going uphill again
				curr_peak = v[i];
			}
		}
	
	}
	# check last point, if still going uphill include it as a maximum
	if (ascending) { pos_peak=cbind(pos_peak,i);n_peaks=n_peaks+1 }
	
	n_peaks
	pos_peak
	length(pos_peak)
	pos_peak=pos_peak[-1]
	#pos_peak
	#length(pos_peak)
	#v[1:30]

	#-----------------------------------------------------
	# remove peaks where v[i]>vmax
	#pos_peak
	#length(pos_peak)
	v[pos_peak]
	if (max(v[pos_peak])>vmax){
		peaks_to_remove<-pos_peak[which(v[pos_peak]>vmax)];peaks_to_remove
		length(peaks_to_remove)
		dset2<-dset2[peaks_to_remove*-1,]
		dset2[1:12,]
		}
	maxi<-max(dset2$v);#maxi

}
#n_int
#length(dset2$lat)
	
# Find removed rows
not_rem<-match(dset2$row_id,dset$row_id);#not_rem
removed<-dset[not_rem*-1,];#removed
vmask<-sda_filter<-character(length(dset$lat));#vmask
vmask[not_rem*-1]<-"removed";#vmask
vmask[not_rem]<-"not";#vmask
extremes<-c(1,2,length(dset$lat)-1,length(dset$lat));#extremes
vmask[extremes]<-"end_location"
vmask
}

