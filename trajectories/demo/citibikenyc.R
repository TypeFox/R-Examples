# https://s3.amazonaws.com/tripdata/201307-citibike-tripdata.zip
# x = read.csv(unz("201307-citibike-tripdata.zip", "2013-07 - Citi Bike trip data.csv"))

library(sp)
library(spacetime)
library(trajectories)

readCBN = function(year = 2013, month = 7) {
	url = "https://s3.amazonaws.com/tripdata/"
	f = "-citibike-tripdata.zip"
	tmp = "CBNTMPxxx.zip"
	res = NULL
	for (y in year) {
		for (m in month) {
			if (m < 10)
				m = paste("0", m, sep = "")
			fname = paste(url, year, m, f, sep = "")
			download.file(fname, tmp, "wget")
			csvname = paste(y, "-", m, " - Citi Bike trip data.csv", sep = "")
			x = read.csv(unz(tmp, csvname))
			x$starttime = as.POSIXct(x$starttime)
			x$stoptime = as.POSIXct(x$stoptime)
			x$birth.year = as.numeric(as.character(x$birth.year))
			x$gender = factor(x$gender, labels = c("unknown", "male", "female"))
			if (is.null(res))
				res = x
			else
				res = rbind(res, x)
		}
	}
	res
}

toTrackLst = function(x) {
	lapply(1:nrow(x), function(i) {
		if (i %% 100 == 0) 
			print(i)
		start = cbind(x[i,"start.station.longitude"], 
			x[i,"start.station.latitude"])
		end = cbind(x[i,"end.station.longitude"], 
			x[i,"end.station.latitude"])
		pt = SpatialPoints(rbind(start, end), 
			CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
		Track(STIDF(pt, c(x[i,"starttime"], x[i,"stoptime"]), 
			data.frame(id=c(x[i,"start.station.id"], x[i,"end.station.id"]))))
	})
}

x = readCBN()
trj = Tracks(toTrackLst(x[1:10,]))
