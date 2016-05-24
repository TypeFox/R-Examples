data.read.phase <- function(path, filename){
	phase <- read.delim(file=paste(path,"/",filename,sep=""))
	phase.tmp <- data.frame(DWD_STAT_ID=phase$DWD_STAT_ID, 
				STAT_NAME=phase$STAT_NAME, 
				STAT_LON=phase$STAT_LON, 
				STAT_LAT=phase$STAT_LAT, 
				STAT_ALT=phase$STAT_ALT, 
				BEGIN_OBS=phase$BEGIN_OBS, 
				END_OBS=phase$END_OBS, 
				REGION=phase$NATURRAUMGRUPPEN_ID,
				PHASE_ID=phase$PHASE_ID, 
				OBS_DAY=phase$OBS_DAY, 
				OBS_YEAR=phase$OBS_YEAR, 
				CHECKED=phase$CHECKED, 
				outlier=phase$outlier)

	phase <- phase.tmp
	rm(phase.tmp)
	return(phase)
}
