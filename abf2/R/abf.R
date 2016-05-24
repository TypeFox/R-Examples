BLOCK_SIZE <- 512
STRING_JUNK_PREFIX <- 44

# read an Axon ABF2 file
# for the moment only gap-free files are supported

abfload <- function ( filename=NULL )
{
	if ( is.null(filename) )
	{
		filename <- file.choose()
	}

	result <- list()
	fp <- file(filename, open="rb")
	result$signature <- readChar(fp, nchars=4, useBytes=TRUE)
	result$filename <- filename
	
	if ( result$signature != "ABF2" )
	{
		close(fp)
		stop("Only ABF2 files supported at present")
	}
	
	class(result) <- "abf2"
	
	vb <- readBin(fp, what="integer", n=4, size=1)
	result$version <- vb[4] + vb[3] * 0.1 + vb[2] * 0.01 + vb[1] * 0.001
	
	result <- defRead(fp, header2Def, result)
	
	result$startTime <- result$startTime * 0.001
	
	# read the section index
	result$sections <- list()
	for ( sect.name in abfSections )
	{
		result$sections[[sect.name]] <- defRead ( fp, sectionInfo, repos=FALSE )
	}
	
	# read the strings section
	seek(fp, where=result$sections$StringSection$blockIndex * BLOCK_SIZE + STRING_JUNK_PREFIX)
	result$strings <- c(readBin(fp, "character", n=result$sections$StringSection$numEntries))
	
	# read the ADC defs
	result$ADC <- list()
	for ( ii in 1:result$sections$ADCSection$numEntries )
	{
		seek(fp, where=result$sections$ADCSection$blockIndex * BLOCK_SIZE + result$sections$ADCSection$bytes * (ii-1))
		result$ADC[[ii]] <- defRead(fp, abfADCInfo, repos=FALSE)
		result$ADC[[ii]]$name <- result$strings[result$ADC[[ii]]$ADCChannelNameIndex]
		result$ADC[[ii]]$units <- result$strings[result$ADC[[ii]]$ADCUnitsIndex]
	}

	# read the protocol block
	seek(fp, where=result$section$ProtocolSection$blockIndex * BLOCK_SIZE)
	result$protocol <- defRead(fp, abfProtocolInfo, repos=FALSE)
	
	# and the tags, if any
	result$rawTags <- list()
	vLevel <- NA
	vUnits <- "" 
	result$tags <- data.frame(time=0, level=vLevel, units=vUnits, comment="", vChange=TRUE)
	if ( result$sections$TagSection$numEntries > 0 )
	{
		for ( ii in 1:result$sections$TagSection$numEntries )
		{
			seek(fp, where=result$sections$TagSection$blockIndex *
					       BLOCK_SIZE + result$sections$TagSection$bytes * (ii-1))
			result$rawTags[[ii]] <- defRead(fp, abfTagInfo, repos=FALSE)
			
			# this is purely empirical based on what I've seen in our own files
			# I have no idea whether there's any more official scheme for this stuff
			rmx <- regexec("[^']*'([^']+)' => ([^ ]+) ([^ ]+)", result$rawTags[[ii]]$comment)
			result$rawTags[[ii]]$is.cmd <- rmx[[1]][1] != -1
			
			tagTime <- result$rawTags[[ii]]$tagTime *
						result$protocol$ADCSequenceInterval *
						1e-6 / length(result$ADC)
			
			if ( result$rawTags[[ii]]$is.cmd )
			{
				mxs <- regmatches(result$rawTags[[ii]]$comment, rmx)[[1]]
				result$rawTags[[ii]]$cmd <- mxs[2]
				result$rawTags[[ii]]$level <- as.numeric(mxs[3])
				result$rawTags[[ii]]$units <- mxs[4]
				
				vLevel <- result$rawTags[[ii]]$level
				vUnits <- result$rawTags[[ii]]$units
				
				result$tags <- rbind(result$tags,
									 data.frame(time=tagTime,
										        level=vLevel,
										        units=vUnits,
										        comment="",
										        vChange=TRUE,
										        stringsAsFactors=FALSE))
			}
			else
			{
				result$comments <- rbind(result$comments,
										 data.frame(time=tagTime,
										            comment=paste(unlist(strsplit(result$rawTags[[ii]]$comment,
										            							  " +")),
										            			  collapse=" "),
										            stringsAsFactors=FALSE))
				result$tags <- rbind(result$tags,
									 data.frame(time=tagTime,
										        level=vLevel,
										        units=vUnits,
										        comment=paste(unlist(strsplit(result$rawTags[[ii]]$comment,
										            						  " +")),
										            		  collapse=" "),
										        vChange=FALSE,
										        stringsAsFactors=FALSE))
			}
		}
	}
	
	# OK, now we need to read the data
    if ( result$dataFormat == 0 )
    {
        dataSz <- 2
        dataType <- "integer"
    }
    else if ( result$dataFormat == 1 )
    {
        dataSz <- 4
        dataType <- "numeric"
    }
    else
    {
        warning(paste("unknown data format", result$dataFormat, "- data section not read"))
        close(fp)
        invisible(result)
    }

    headOffset <- result$sections$DataSection$blockIndex * BLOCK_SIZE
    #si <- result$protocol$ADCSequenceInterval

    if ( result$protocol$operationMode == 3 )
    {
        # gap free mode - this is our primary goal
        #result$dataPtsPerChannel <- result$sections$DataSection$numEntries / length(result$ADC)
        seek(fp, where=headOffset)
        result$rawdata <- readBin(fp, what=dataType, size=dataSz, n=result$sections$DataSection$numEntries)
        result$traces <- matrix(result$rawdata, nrow=length(result$ADC))
        result$s <- 0:(dim(result$traces)[2]-1) * (result$protocol$ADCSequenceInterval * 1e-6)
    
        # scale int data by ADC settings
        if ( result$dataFormat == 0 )
        {
            for ( ii in 1:length(result$ADC) )
            {
                adc <- result$ADC[[ii]]
                result$ADC[[ii]]$scaleFactor <- result$protocol$ADCRange /
                                                (result$protocol$ADCResolution *
                                                adc$instScaleFactor *
                                                adc$signalGain *
                                                adc$ADCProgGain *
                                                adc$teleAddGain)
                                                                
                result$ADC[[ii]]$offset <- adc$instOffset - adc$signalOffset
            
                result$traces[ii,] <- result$traces[ii,] * result$ADC[[ii]]$scaleFactor + result$ADC[[ii]]$offset
            }
        }
    }
    else
    {
        warning(paste("unsupported data mode", result$protocol$operationMode, "- data section not read"))
    }
	
	close(fp)
	invisible(result)
}

defRead <- function ( fp, def, result=NULL, origin=0, repos=TRUE )
{
	if ( is.null(result) )
		result <- list()
	
	for ( ii in 1:length(def$field) )
	{
		if (repos) { seek(fp, where=def$offset[ii]+origin) }
		if ( def$type[ii]=="uint32" )
		{
			# R does not support unsigned long, so we have to fudge this to double
			lo <- readBin(fp, what="integer", size=2, signed=FALSE)
			hi <- readBin(fp, what="integer", size=2, signed=FALSE)
			result[[def$field[ii]]] <- 65536.0 * hi + lo
		}
		else if ( def$type[ii]=="string" )
		{
			result[[def$field[ii]]] <- readChar(fp, nchars=def$bytes[ii], useBytes=TRUE)
		}
		else if ( def$type[ii]=="skip" )
		{
			seek(fp, where=def$bytes[ii], origin="current")
		}
		else
		{
			result[[def$field[ii]]] <- readBin(fp, what=def$type[ii], size=def$bytes[ii])
		}
	}
	
	invisible(result)
}

plot.abf2 <- function ( x,
						adc=1,
						time=NULL,
						pts=1000,
						type="s",
						col=1,
						vtags=TRUE,
						ctags=TRUE,
						vcol="grey",
						ccol="pink",
						vlty="dotted",
						clty="dashed",
						xlab="Time (s)",
						ylab=NULL,
						vText="bottom",
						vtunits=FALSE,
						vtcol=4,
						vtpad=0.1,
						... )
{
	if ( length(time) == 0 )
	{
		time <- c(0,max(x$s))
	}
	else if ( length(time) == 1 )
	{
		time <- c(0,time)
	}
	
	dt <- x$s[2] - x$s[1]
	nsamps <- ceiling((time[2] - time[1])/dt)
	
	if ( length(pts)==0 || pts > nsamps )
	{
		idx <- which(x$s >= time[1] & x$s <= time[2])
	}
	else
	{
		idx <- floor((time[1]/dt) + (1:pts * nsamps/(pts-1)))
		idx <- idx[idx <= length(x$s)]
	}
	
	if ( length(ylab) == 0 )
	{
		ylab=paste(x$ADC[[adc]]$name, " (", x$ADC[[adc]]$units, ")", sep="")
	}
	
	lo <- min(x$traces[adc,idx])
	hi <- max(x$traces[adc,idx])
	
	if ( vtags && sum(x$tags$vChange) > 0 )
	{
		vt <- x$tags[x$tags$vChange & x$tags$time >= time[1] & x$tags$time <= time[2],]
		
		if ( length(vt$level) > 0 && length(vText) > 0 )
		{
			pad <- (hi - lo) * vtpad
			if ( vText=="bottom" )
			{
				lo <- lo - pad
				vty <- rep(c(lo + pad/4, lo + 3 * pad/4), length.out=length(vt$level))
			}
			else
			{
				vty <- rep(c(hi + pad/4, hi + 3 * pad/4), length.out=length(vt$level))
				hi <- hi + pad
			}
		}
	}
	else
	{
		vt <- NULL
	}
	
	if ( ctags && sum(! x$tags$vChange) > 0 )
	{
		ct <- x$tags[!x$tags$vChange & x$tags$time >= time[1] & x$tags$time <= time[2],]
	}
	else
	{
		ct <- NULL
	}
		
	plot( x=x$s[idx], y=x$traces[adc,idx], type=type, col=col, xlab=xlab, ylab=ylab, ylim=c(lo,hi), ... )
	
	if ( length(vt$level > 0 ) )
	{
		abline(v=vt$time, col=vcol, lty=vlty)
		
		if ( length(vText) > 0 )
		{
			if ( vtunits )
			{
				labels <- paste(vt$level, vt$units)
			}
			else
			{
				labels <- vt$level
			}
			
			text(x=vt$time, y=vty, labels=labels, col=vtcol)
		}
	}
	
	if ( length(ct$time) > 0 )
	{	
		abline(v=ct$time, col=ccol, lty=clty)
	}
}

# definition for ABF2 header fields we're interested in
# note that the signature and version number are not included
# since we obtain that info beforehand
header2Def <- data.frame(
    field=c("headerSize", "episodes", "startDate", "startTime",
    		"stopwatchTime", "fileType", "dataFormat", "nScan",
    		"CRCEnable", "fileCRC", "fileGUID", "creatorVersion",
    		"creatorNameIndex", "modifierVersion", "modifierNameIndex", "protocolPathIndex"),
    offset=c(8, 12, 16, 20,
    		 24, 28, 30, 32,
    		 34, 36, 40, 56,
    		 60, 64, 68, 72),
    type=c("uint32", "uint32", "uint32", "uint32",
           "uint32", "integer", "integer", "integer",
           "integer", "uint32", "uint32", "uint32",
           "uint32", "uint32", "uint32", "uint32"),
    bytes=c(4, 4, 4, 4,
    		4, 2, 2, 2,
    		2, 4, 4, 4,
    		4, 4, 4, 4),
    stringsAsFactors=FALSE
)

# sections
abfSections <- c(
	"ProtocolSection",
	"ADCSection",
	"DACSection",
	"EpochSection",
	"ADCPerDACSection",
	"EpochPerDACSection",
	"UserListSection",
	"StatsRegionSection",
	"MathSection",
	"StringSection",
	"DataSection",
	"TagSection",
	"ScopeSection",
	"DeltaSection",
	"VoiceTagSection",
	"SynchArraySection",
	"AnnotationSection",
	"StatsSection"
)

sectionInfo <- data.frame(
	field=c("blockIndex", "bytes", "numEntries"),
	type=c("uint32", "uint32", "integer"),
	bytes=c(4,4,8),
	stringsAsFactors=FALSE
)

abfProtocolInfo <- data.frame(
	field=c("operationMode", "ADCSequenceInterval", "fileCompressionEnabled", "unused",
			"fileCompressionRatio", "synchTimeUnit", "secondsPerRun", "samplesPerEpisode",
			"pretriggerSamples", "episodesPerRun", "runsPerTrial", "nTrials",
			"averagingMode", "undoRunCount", "firstEpInRun", "triggerThreshold",
			"triggerSource", "triggerAction", "triggerPolarity", "scopeOutputInterval",
			"episodeStartToStart", "runStartToStart", "averageCount", "trialStartToStart",
			"autoTriggerStrategy", "firstRunDelay", "channelStatsStrategy", "samplesPerTrace",
			"startDisplayNum", "finishDisplayNum", "showPNRawData", "statsPeriod",
			"statsMeasurements", "statsSaveStrategy", "ADCRange", "DACRange",
			"ADCResolution", "DACResolution", "experimentType", "manualInfoStrategy",
			"commentsEnable", "fileCommentIndex", "autoAnalyseEnable", "signalType",
			"digitalEnable", "activeDACChannel", "digitalHolding", "digitalInterEpisode",
			"digitalDACChannel", "digitalTrainActiveLogic", "statsEnable", "statsClearStrategy",
			"levelHysteresis", "timeHysteresis", "allowExternalTags", "averageAlgorithm",
			"averageWeighting", "undoPromptStrategy", "trialTriggerSource", "statsDisplayStrategy",
			"externalTagType", "scopeTriggerOut", "LTPType", "alternateDACOutputState",
			"alternateDigitalOutputState", "cellID1", "cellID2", "cellID3",
			"digitizerADCs", "digitizerDACs", "digitizerTotalDigOuts", "digitizerSynchDigOuts",
			"digitizerType"),
	type=c("integer", "numeric", "logical", "skip",
			"uint32", "numeric", "numeric", "integer",
			"integer", "integer", "integer", "integer",
			"integer", "integer", "integer", "numeric",
			"integer", "integer", "integer", "numeric",
			"numeric", "numeric", "integer", "numeric",
			"integer", "numeric", "integer", "integer",
			"integer", "integer", "integer", "numeric",
			"integer", "integer", "numeric", "numeric",
			"integer", "integer", "integer", "integer",
			"integer", "integer", "integer", "integer",
			"integer", "integer", "integer", "integer",
			"integer", "integer", "integer", "integer",
			"integer", "integer", "integer", "integer",
			"numeric", "integer", "integer", "integer",
			"integer", "integer", "integer", "integer",
			"integer", "numeric", "numeric", "numeric",
			"integer", "integer", "integer", "integer",
			"integer"),
	bytes=c(2, 4, 1, 3,
			4, 4, 4, 4,
			4, 4, 4, 4,
			2, 2, 2, 4,
			2, 2, 2, 4,
			4, 4, 4, 4,
			2, 4, 2, 4,
			4, 4, 2, 4,
			4, 2, 4, 4,
			4, 4, 2, 2,
			2, 4, 2, 2,
			2, 2, 2, 2,
			2, 2, 2, 2,
			2, 4, 2, 2,
			4, 2, 2, 2,
			2, 2, 2, 2,
			2, 4, 4, 4,
			2, 2, 2, 2,
			2),
			stringsAsFactors=FALSE
)

abfMathInfo <- data.frame(
	field=c("mathEnable", "mathExpression", "mathOperatorIndex", "mathUnitsIndex",
			"mathUpperLimit", "mathLowerLimit", "mathADCNum1", "mathADCNum2",
			"unused", "mathK1", "mathK2", "mathK3",
			"mathK4", "mathK5", "mathK6"),
	type=c("integer", "integer", "uint32", "uint32",
			"numeric", "numeric", "integer", "integer",
			"skip", "numeric", "numeric", "numeric",
			"numeric", "numeric", "numeric"),
	bytes=c(2, 2, 4, 4,
			4, 4, 2, 2,
			16, 4, 4, 4,
			4, 4, 4),
	stringsAsFactors=FALSE
)

abfADCInfo <- data.frame(
	field=c("ADCNum", "teleEnable", "teleInstrument", "teleAddGain",
			"teleFilter", "teleMembraneCap", "teleMode", "teleAccResist",
			"ADCPtoLChannelMap", "ADCSamplingSeq", "ADCProgGain", "ADCDispAmp",
			"ADCDispOffset", "instScaleFactor", "instOffset", "signalGain",
			"signalOffset", "signalLowpass", "signalHighpass", "lowpassType",
			"highpassType", "postprocLowpass", "postprocLowpassType", "enabledDuringPN",
			"statsChannelPolarity", "ADCChannelNameIndex", "ADCUnitsIndex"),
	type=c("integer", "integer", "integer", "numeric",
			"numeric", "numeric", "integer", "numeric",
			"integer", "integer", "numeric", "numeric",
			"numeric", "numeric", "numeric", "numeric",
			"numeric", "numeric", "numeric", "integer",
			"integer", "numeric", "integer", "logical",
			"integer", "integer", "integer"),
	bytes=c(2, 2, 2, 4,
			4, 4, 2, 4,
			2, 2, 4, 4,
			4, 4, 4, 4,
			4, 4, 4, 1,
			1, 4, 1, 1,
			2, 4, 4),
	stringsAsFactors=FALSE
)

abfTagInfo <- data.frame(
	field=c("tagTime", "comment", "tagType", "annotIndex"),
	type=c("integer", "string", "integer", "integer"),
	bytes=c(4, 56, 2, 2),
	stringsAsFactors=FALSE
)

split.abf2 <- function ( x, f=NULL, drop=FALSE, adc=1, lag=0.3, ... )
{
	if ( class(x) != "abf2" )
	{
		stop("can only split traces for class 'abf2'")
	}
	
	if ( adc > length(x$ADC) )
	{
		stop(paste("unknown ADC channel:", adc))
	}
	
	if ( length(x$tags$time) < 2  )
	{
		stop("no tags to split")
	}
	
	result <- list()
	for ( ii in 1:(length(x$tags$time)) )
	{
		hi <- ifelse(ii == length(x$tags$time), max(x$s), x$tags$time[ii+1])
		ind <- which((x$s >= x$tags$time[ii] + lag) & (x$s <= hi))
		result[[ii]] <- list()
		result[[ii]]$trace <- x$traces[adc,ind]
		result[[ii]]$s <- x$s[ind]
		result[[ii]]$level <-x$tags$level[ii]
		result[[ii]]$index <- ii
		result[[ii]]$comment <- x$tags$comment[ii]
		class(result[[ii]]) <- "abf2split"
	}
	
	invisible(result)
}

print.abf2 <- function ( x, ... )
{
	summary.abf2(x, tags=FALSE, ... )
}

summary.abf2 <- function ( object, tags=TRUE, ... )
{
	if ( class(object) != "abf2" )
	{
		stop("object of class 'abf2' required")
	}
	
	cat(paste("\n", ifelse(object$dataFormat,"Floating point ", "Integer "), "ABF2 loaded from ", object$filename, "\n", sep=""))
	
	dmy <- object$startDate
	y <- floor(dmy/1e4)
	dm <- dmy - y * 1e4
	mo <- floor(dm/1e2)
	d <- dm - mo * 1e2
	
	hms <- object$startTime
	h <- floor(hms/3600)
	ms <- hms - h * 3600
	mi <- floor(ms/60)
	s <- ms - mi * 60
	
	dur <- object$s[length(object$s)]
	durh <- floor(dur/3600)
	durms <- dur - durh * 3600
	durm <- floor(durms/60)
	durs <- durms - durm * 60
	
	cat(paste("Recorded ", y, "-",
			  ifelse(mo<10,"0",""), mo, "-",
			  ifelse(d<10,"0",""), d, " ",
			  hms2str(h, mi, s), "\n",
			  sep=""))
	cat(paste("Duration ", hms2str(durh, durm, durs),
			  " (", formatC(dur, digits=2, format="f"), " s)\n\n",
			  sep=""))
			  
	intvl <- object$s[2] - object$s[1]
	hz <- floor(1/intvl)
	nsamp <- dim(object$traces)[2]
	cat(paste(length(object$ADC), " channel", ifelse(length(object$ADC) > 1, "s", ""), ", sample interval ", intvl * 1e3, " ms (", hz, " Hz), ", nsamp, " samples/channel\n\n",
			  sep=""))
	
	for ( ii in 1:length(object$ADC) )
	{
		cat( paste( " [", ii, "] #", object$ADC[[ii]]$ADCNum, ": ",
			        object$ADC[[ii]]$name, " (", object$ADC[[ii]]$units, "; ",
			        object$ADC[[ii]]$teleFilter, " Hz tele filter, x", object$ADC[[ii]]$teleAddGain, " tele gain)\n",
			        sep="") )
	}
	
	cat("\n")
	cat(paste(length(object$rawTags), " tags, ", sum(object$tags$vChange) - 1, " detected voltage changes\n", sep=""))
	
	if ( tags && length(object$rawTags) > 0 )
	{
		cat("\n")
		tt <- NULL
		for ( ii in 1:length(object$rawTags) )
		{
			tt <- c(tt, object$rawTags[[ii]]$tagTime * object$protocol$ADCSequenceInterval * 1e-6 / length(object$ADC))
		}
		
		tts <- formatC(tt, digits=2, format="f")
		pad <- max(nchar(tts)) - nchar(tts) + 1
		
		for ( ii in 1:length(object$rawTags) )
		{
			pre <- 1 + floor(log(length(object$rawTags), 10)) - floor(log(ii,10))
			cat(paste(rep(" ", pre), collapse=""))
			cat(paste("[", ii, "]", sep=""))
			cat(paste(rep(" ", pad[ii]), collapse=""))
			cat(tts[ii])
			cat(paste(": ", object$rawTags[[ii]]$comment, "\n", sep=""))
		}
	}
	
	cat("\n")
}

hms2str <- function(h,m,s)
{
	return( paste( ifelse(h < 10, "0", ""), floor(h), ":",
				   ifelse(m < 10, "0", ""), floor(m), ":",
				   ifelse(s < 10, "0", ""), formatC(s, digits=2, format="f"),
				   sep=""))
}

plot.abf2split <- function ( x,
							 pts=1000,
							 time=NULL,
							 type="s",
							 col=1,
							 xlab="Time (s)",
							 ylab="Current (pA)",
							 main=paste("Segment ", x$index, " at ", x$level, " mV  [", x$comment, "]", sep="" ),
							 ... )
{
	if ( length(time) == 0 )
	{
		time <- c(min(x$s),max(x$s))
	}
	else if ( length(time) == 1 )
	{
		time <- c(min(x$s),time)
	}
	
	dt <- x$s[2] - x$s[1]
	nsamps <- ceiling((time[2] - time[1])/dt)
	
	if ( length(pts)==0 || pts > nsamps )
	{
		idx <- which(x$s >= time[1] & x$s <= time[2])
	}
	else
	{
		idx <- floor((time[1] - x$s[1])/dt) + 1:pts * floor(nsamps/(pts-1))
		idx <- idx[idx <= length(x$s)]
	}

	plot( x=x$s[idx], y=x$trace[idx], type=type, col=col, xlab=xlab, ylab=ylab, main=main, ... )	
}

# plot multiple segments from the same trace
multiplot <- function ( x,
						adc=1,
						duration=1,
						start=0,
						pts=1000,
						type="s",
						single.col=1,
						local.col=TRUE,
						gutter=NULL,
						gutter.prop=0.1,
						labels=NULL,
						rotate.labels=90,
						time.scale=0.2,
						time.scale.label=paste(time.scale*1000,"ms",sep=""),
						trace.scale=5,
						trace.scale.label=paste(trace.scale,"pA",sep=""),
						scale.col="grey50",
						xinset=NULL,
						xinset.prop=1/20,
						... )
{
	if ( class(x) != "abf2" )
	{
		stop("argument must be of class 'abf2'")
	}
	
	n <- length(start)
	dt <- x$s[2] - x$s[1]
	nsamps <- ceiling(duration/dt)
	
	if ( length(pts)==0 || pts > nsamps )
	{
		idx <- 1:nsamps
	}
	else
	{
		idx <- 1:pts * floor(nsamps/(pts-1))
	}
	
	t <- 1:length(idx) * duration/length(idx)
	
	traces <- matrix(0, nrow=n, ncol=length(idx))
	spans <- vector(mode="numeric", length=n)
	offsets <- rep(0, n)
	voltages <- vector(mode="numeric", length=n)
	
	# gather trace details	
	for ( ii in 1:n )
	{
		traces[ii,] <- x$traces[adc, idx + floor(start[ii]/dt)]
		traces[ii,] <- traces[ii,] - min(traces[ii,])
		spans[ii] <- max(traces[ii,])
		
		if ( !is.null(x$tags[x$tags$vChange]) )
		{
			vix <- which(x$tags[x$tags$vChange]$time <= start[ii])
			if ( length(vix) > 0 )
			{
				voltages[ii] <- x$tags[x$tags$vChange]$level[vix[length(vix)]]
			}
			else
			{
				voltages[ii] <- NA
			}
		}
	}
	
	if ( is.null(gutter) )
	{
		gutter <- gutter.prop * max(spans)
	}
	
	if ( is.null(xinset) )
	{
		xinset <- duration * xinset.prop
	}
	
	ylim <- c(0, sum(spans) + (n-1) * gutter)
	xlim <- c(-xinset, duration)
	
	if ( local.col )
	{
		clut <- rainbow(n, end=5/6)
	}
	else
	{
		clut <- rep(single.col, n)
	}
	
	if ( is.null(labels) )
	{
		labels <- voltages
	}
	else
	{
		labels <- rep(labels, length.out=n)
	}
	
	plot( x=t,
		  y=traces[1,],
		  type=type,
		  col=clut[1],
		  ylab="",
		  xlab="",
		  ylim=ylim,
		  xlim=xlim,
		  xaxt="n",
		  yaxt="n",
		  bty="n",
		  ... )
	text( -xinset,
		  spans[1]/2,
		  as.character(labels[1]),
		  srt=rotate.labels )
	
	if ( n > 1 )
	{
		for ( ii in 2:n )
		{
			offsets[ii] <- offsets[ii-1] + spans[ii-1] + gutter

			lines( x=t,
				   y=traces[ii,] + offsets[ii],
				   col=clut[ii],
				   type=type )
			text( -xinset,
				  offsets[ii]+spans[ii]/2,
				  as.character(labels[ii]),
				  srt=rotate.labels )
		}
	}
	
	x.o <- -2 * xinset
	y.o <- -0.1 * ylim[2]
		
	old.par <- par(xpd=TRUE,lwd=4)
		
	lines(c( x.o, x.o, x.o + time.scale ), c( y.o + trace.scale, y.o, y.o), col=scale.col)
	text( x.o + (time.scale/2), y.o * 1.4, time.scale.label, col=scale.col )
	text( x.o * 1.3, y.o + (trace.scale/2), trace.scale.label, srt=rotate.labels, col=scale.col )
	par(old.par)
}
