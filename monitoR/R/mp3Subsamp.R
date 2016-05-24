# Function to subsample recordings for shorter surveys, e.g. 10 minutes from the top of every hour
# Requires mp3splt, libmp3splt to be installed separately
# Written by J. Katz 08-May 2013
# Modified: 2014 MAR 29

# To install mp3splt on windows you need to download the zip or exe and then add the path to your exe.
# Search for cmd, then type this into the cmd window--replacing the path below with your actual path.
# path C:\Program Files\mp3splt;%PATH%
# 

mp3Subsamp <- function(
    files,              # Vector or list of files to transfer, or specify 'from'
    from='.',           # Entire SD card or other directory with field recordings, or specify 'files'
    to,                 # Destination directory
    csv.dir=to,         # Directory where table should be saved
    csv.name,           # Name of csv file
    duration=600,       # Duration (seconds) of surveys to pull
    mins.between=50,    # Number of minutes between surveys
    index='hour',       # c('hour', 'time0') Work from top of every hour or time 0 of the file
    loc.prefix,         # 6-character location prefix, e.g MABI01
    CardRecorderID=NA,  # Associates surveys with a recorder and SD card
    kbps=128,           # c(128, 192, 256, 320) mp3 bitrate
    samp.rate=44100,    # mp3 sample rate
    channels=2,         # Stereo & JntStereo=2, Mono=1
    split=TRUE          # Set to FALSE for metadata only
    ){
    
    if (length(system(command="mp3splt -h", intern=TRUE))<10) stop("'mp3splt' not found; download mp3splt at\nhttp://mp3splt.sourceforge.net/mp3splt_page/home.php\n")

    if(missing(csv.name)) csv.name <- paste(loc.prefix, '_', as.character(Sys.time(), format='%Y-%m-%d'), '.csv', sep='')

    if (nchar(loc.prefix) != 6) {stop(paste('loc.prefix must be 6 characters, got', loc.prefix))}

    # Create "to" directory if it doesn't exist
    if(missing(to)) stop('No destination directory specified')
    if(!file.exists(to)) dir.create(to)

    if(missing(files)){ # if moving a whole directory:
        files <- list.files(from, full.names=TRUE, pattern='[mMpP3]')
        survey.names <- list.files(from, pattern='[mMpP3]')
    } else { # if only moving a few files:
        survey.full.names <- strsplit(as.character(files), "/")
        survey.names <- unlist(lapply(survey.full.names, function(x) x[length(x)]))
    }
    surveymeta.l <- lapply(as.list(files), function(x) mp3Subsamp.one(recording=x, duration=duration, mins.between=mins.between, index=index, loc.prefix=loc.prefix, to=to, CardRecorderID=CardRecorderID, kbps=kbps, samp.rate=samp.rate, channels=channels, split=split))
    surveymeta.dat <- rbindf(surveymeta.l)
    write.csv(surveymeta.dat, paste0(csv.dir, "/", csv.name))
    invisible(surveymeta.dat)
    }

mp3Subsamp.one <- function(
    recording,          # Name of recording from which to draw surveys
    duration,           # Duration (seconds) of surveys to pull
    mins.between,       # Number of minutes between surveys
    index,              # c('hour', 'time0') Work from top of every hour or time 0 of the file
    loc.prefix,         # 6-character location prefix, e.g MABI01
    to,                 # Destination directory for new subsurveys
    CardRecorderID,     # Associates surveys with a recorder and SD card
    kbps,               # mp3 bitrate
    samp.rate,          # mp3 sample rate
    channels,           # Stereo & JntStereo=2, Mono=1
    split               # Set to FALSE for metadata only
    ){

    # Get mdate
    file.mdate <- file.info(recording)$mtime
    bps <- kbps*1000/8
    size <- file.info(recording)$size
    rec.dur <- floor(size/bps)
    # Determine file start time
    stime <- file.mdate-rec.dur	
    # Determine time in file of first survey; if index == 'hour' it is the first hour, otherwise it is the first segment of recording
    if(tolower(index) %in% c('hour', 'h', 'hr')) {wait <- 3600-as.numeric(stime)%%3600
    } else if(tolower(index) %in% c('time0', 't', 'time', '0')) wait <- 0
    # Make from and to lists for mp3splt
    surv.start <- seq(wait, rec.dur, by=duration+mins.between*60)
    surv.end <- surv.start+duration
    # predict the new mdates and file names for each subsurvey
    time.remain <- rec.dur-surv.end
    new.file.mdates <- as.character(file.mdate-time.remain, format='%Y-%m-%d %H%M%S %Z')
    new.mdates <- as.character(file.mdate-time.remain, format='%Y-%m-%d_%H%M%S_%Z')
    filenames <- paste(loc.prefix, '_', new.mdates, sep="")
    # Format seconds to mm.ss.xx, where xx are hundredths, for mp3splt 
    surv.start <- lapply(surv.start, function(x) {
        if(x>=60) { 
            min <- x%/%60
            sec <- floor(x%%60)+round(x%%floor(x), 2)            
            if(nchar(sec)<5 && round(x%%floor(x), 2)>0) sec <- paste0('0', sec)
            from <- as.numeric(paste0(min, '.', sec))
        } else if (x<10) {
            from <- as.numeric(paste0('00.0', round(x, 2)))
        } else from <- as.numeric(paste0('00.', round(x, 2)))
        return(from)
        })
    surv.end <- lapply(surv.end, function(x) {
        if (x>=60) {
            min <- x%/%60
            sec <- floor(x%%60)+round(x%%floor(x), 2)
            if(nchar(sec)<5 && round(x%%floor(x), 2)>0) sec <- paste0('0', sec)
            to <- as.numeric(paste0(min, '.', sec))
            if(to*60>rec.dur) to <- 'EOF'
        } else if (to<10) {
            to <- as.numeric(paste0('00.0', round(x, 2)))
        } else to <- as.numeric(paste0('00.', round(x, 2)))
        return(to)
        })
    # Call mp3splt to cut up the mp3 file
    if(split){
        if(missing(to)) to <- getwd()
        num.survs <- as.list(1:length(surv.start))
        lapply(num.survs, function(x) system(paste0("mp3splt -Q -d ", to, " ", recording, " ", surv.start[[x]], " ", surv.end[[x]], " -o ", filenames[x])) )
    #    lapply(num.survs, function(x) system(paste0("mp3splt -Q -d '", to, "' '", recording, "' ", surv.start[[x]], " ", surv.end[[x]], " -o '", filenames[x], "'")) ) # Quotes don't work in windows?
    }
    filenames <- paste0(filenames, '.mp3')

    return(data.frame(fldOriginalDateModified=new.file.mdates, fldOriginalRecordingName=recording, fldSurveyName=filenames, fldRecordingFormat='mp3', fkCardRecorderID=CardRecorderID, fldSurveyLength=duration, fldSampleRate=samp.rate, fldBitsperSample=16, fldChannels=channels))
}





    
