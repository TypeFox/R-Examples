flux.calib <-
function(dat, columns, calib, format = "%Y-%m-%d %H:%M:%S", window = 3, buffer = 1000, n.cg = 4, rl.backup = 20, attach = FALSE){
	# defining the function which does the work
	flux.cal <-
	function(conz.dat, calib, format = "%Y-%m-%d %H:%M:%S", window=3, buffer = 1000, n.cg = 4){
		# do not allow NA measurements in conz.dat
		conz.dat <- conz.dat[!is.na(conz.dat[,2]),]
		# extract date
		m.date <- strptime(conz.dat[1,1], format = format)
		dts.cal <- strptime(calib[,1], format = format)
		na.omit <- !is.na(dts.cal)
		dts.cal <- dts.cal[na.omit]
		calib <- calib[na.omit,]
		# extract calibration gas measurements according to the date of the measurement
		# of the ghg and a window width window (hours) around it
		# first make seconds window because seconds are the primary unit for datetime objects
		window <- window*60*60/2
		calib <- calib[(dts.cal >= (m.date-(window-60))) & (dts.cal <= (m.date+window)), 2]
		# omit 0 concentrations
		calib <- calib[calib>1]		
		# check whether enough concentrations are in calib
		if(length(calib) <= n.cg){
			range.lim <- rl.backup
			warning("Not enough calibration gas measurements", call. = FALSE)
		}
		else{
			# create an index from the grouping of the calibration gas measurements (via clustering)
			cin <- cutree(hclust(dist(calib)), n.cg)
			# omit calibration gas concentrations that are too far away from measured concentrations
			# to provide helpful ranges
			sel <- rowSums(as.matrix(dist(c(range(conz.dat[,2]), calib)))[-c(1:2),1:2] < buffer) > 0
			# calculate range limits (standard deviation of the calibration gas 
			# measurements) per calibration gas  
			range.lims <- as.vector(by(data.frame(calib, sel), cin, function(x) with(x, sd(calib[sel]))))
			# take only the good ones (non that are more than 500 apart)
			#tmp <- as.matrix(dist(c(range(conz.dat[,2]),range.lims)))[-c(1:2),1:2]
			#which <- apply(tmp,1,function(x) sum(x>500)>0)
			# calculate average range limits across all included calibration gases
			range.lim <- mean(range.lims, na.rm=TRUE)
		}	
		return(range.lim)
	}
	# actually do the work
	# extract the needed columns from calib
	calib <- calib[,columns]
	ghg.lim <- sapply(dat$tables, function(x) flux.cal(x[,columns], calib[,columns], format = format, window = window, buffer=buffer, n.cg=n.cg))
	# if NA's result fall back to rl.backup
	ghg.lim[is.na(ghg.lim)] <- ifelse(!is.null(rl.backup), rl.backup, min(ghg.lim, na.rm=TRUE))
	if(attach){
		if(length(grep("CO2", columns))!=0){
			for(i in c(1:length(dat$tables))){
				dat$tables[[i]]$CO2.rl <- ghg.lim[i]
			}
		}
		if(length(grep("CH4", columns))!=0){
			for(i in c(1:length(dat$tables))){
				dat$tables[[i]]$CH4.rl <- ghg.lim[i]
			}
		}
		if(length(grep("N2O", columns))!=0){
			for(i in c(1:length(dat$tables))){
				dat$tables[[i]]$N2O.rl <- ghg.lim[i]
			}
		}					
	}
	else{
		dat <- ghg.lim
	}
	return(dat)
}