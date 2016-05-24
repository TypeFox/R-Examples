############################################################################
# R functions for pawacc 1.2
# Marco Geraci, 20 June 2014
############################################################################

# startup message

".onAttach" <- function(lib, pkg) {
    if(interactive() || getOption("verbose"))
	packageStartupMessage(sprintf("Package %s (%s) loaded. Type citation(\"%s\") on how to cite this package\n", pkg,
		packageDescription(pkg)$Version, pkg))
}

# Generics

markwear <- function(object, value, which = "counts",  rescale.epoch = 60, nz = 0, keep.error = FALSE) UseMethod("markwear")
markpa <- function(object, value, which = "counts", rescale.epoch = 60, labels = NULL, extreme = NULL, keep.error = FALSE) UseMethod("markpa")
markbouts <- function(object, value, which = "counts", bts = c(0,10,20,Inf), rescale.epoch = 60, collapse.by = "%Y-%m-%d", value.labels = NULL, bouts.labels = NULL, extreme = NULL, keep.error = FALSE) UseMethod("markbouts")
collapse <- function(...) UseMethod("collapse")

# Read accelerometer files

readAccDir <- function(path, model, ext = "dat", counts.pos = 1, tz = "Europe/London", sparse = FALSE, fault = 32767, save = TRUE, compress = "gzip", compression_level = 6, ...){

ptm <- proc.time()
 
# read directory content
FILES <- list.files(path, all.files = FALSE, ...)

	if(length(FILES) == 0)
		stop("Directory is empty")

# check names
nChars <- nchar(FILES)
check <- substr(FILES, nChars - 2, nChars) == ext

	if(!all(check))
		warning("Some files do not match specified extension 'ext' and have been ignored")
	if(sum(check) == 1)
		stop("Directory specified contains 1 file. Please see '?gt1mAccFile' instead")
	if(sum(check) == 0)
		stop("Invalid directory")

FILES <- FILES[check]
FILEID <- unlist(strsplit(FILES, paste(".", ext, sep = "")))

# return an info list

accFileList <- list(FILES = FILES, FILEID = FILEID, path = path, model = model, counts.pos = counts.pos, tz = tz, sparse = sparse, fault = fault)

out <- switch(model,
		gt1m = gt1mAccDir(accFileList, save, compress = compress, compression_level = compression_level),
		gt3x = gt3xAccDir(accFileList, save, compress = compress, compression_level = compression_level))

cat("Total elapsed time", (proc.time() - ptm)[3], "seconds\n")
		
return(out)		
		
}

gt1mAccDir <- function(accFileList, save, compress = "gzip", compression_level = 6){

FILES <- accFileList$FILES
FILEID <- accFileList$FILEID
path <- accFileList$path
model <- accFileList$model
N <- length(FILES)
newEnv <- new.env()
tz <- accFileList$tz

# determine saving mode

	if(is.logical(save)){
		if(save){
			saveDir <- paste(path, "/accRfiles", sep = "")
			if(file.exists(saveDir)){
				cat(paste("Directory ", saveDir, " already exists. \n Do you still want to save R files in that folder (y/n)? \n", sep = ""))
				ans <- readLines(stdin(),1)
				if(as.character(ans) == "n")
					stop("Data processing stopped by user")
						else if(as.character(ans) != "y") stop("Ambiguous answer. Process halted.")
				} else {
				dir.create(saveDir)
				cat("Directory ", saveDir, " created \n")
				}
		} else {saveDir <- NULL}
	}

	if(is.character(save)){
		saveDir <- save
		if(file.exists(saveDir)){
			cat(paste("Directory ", saveDir, " already exists. \n Do you still want to save R files in that folder (y/n)? \n", sep = ""))
			ans <- readLines(stdin(),1)
				if(as.character(ans) == "n")
					stop("Data processing stopped by user")
				else if(as.character(ans) != "y") stop("Ambiguous answer. Process halted.")
		}
		if(!file.exists(saveDir)){
			cat("Directory for output does not exist. \n Do you want to create it (y/n)? \n")
			ans <- readLines(stdin(),1)
				if(as.character(ans) == "n")
					stop("Data processing stopped by user")
				else if(as.character(ans) != "y") stop("Ambiguous answer. Process halted.")
					else if(as.character(ans) == "y"){
						dir.create(saveDir)
						cat("Directory ", saveDir, " created \n")
					}
		}
	}

# Loop over files

summaryFile <- data.frame(fileid = "", serial = "", nobs = 0, epoch = 0, mode = NA, ts_start = as.POSIXlt("2000/01/01 00:00:00"), tz = "", voltage = NA, ts_dl = as.POSIXlt("2000/01/01 00:00:00"), stringsAsFactors = FALSE)

pb <- winProgressBar(title = "progress bar", min = 0, max = N, width = 300)
	if(!is.null(saveDir)){
		for(i in 1:N){
			object <- gt1mAccFile(file = FILES[i], path = path, fileid = FILEID[i], counts.pos = accFileList$counts.pos, tz = tz, sparse = accFileList$sparse, fault = accFileList$fault)
			summaryFile[i,] <- object$info
			assign(x = FILEID[i], value = object, envir = newEnv)
			save(list = c(FILEID[i]), file = paste(saveDir, "/", FILEID[i], ".Rdata", sep =""), envir = newEnv, compress = compress, compression_level = compression_level)
			remove(list = c(FILEID[i]), envir = newEnv)
			setWinProgressBar(pb, i, title = paste(round(i/N*100, 0), "% done"))
		}
		close(pb)
		cat("Output is ready \n")
		out <- summaryFile
		class(out) <- c("acclist", "gt1m", "summaryFile")
	} else {
		out <- vector("list", N)
		for(i in 1:N){
			out[[i]] <- gt1mAccFile(file = FILES[i], path = path, fileid = FILEID[i], counts.pos = accFileList$counts.pos, tz = tz, sparse = accFileList$sparse, fault = accFileList$fault)
			summaryFile[i,] <- out[[i]]$info
			setWinProgressBar(pb, i, title = paste(round(i/N*100, 0), "% done"))
		}
		close(pb)
		names(out) <- FILEID
		cat("Output is ready \n")
		attr(out, "info") <- summaryFile
		class(out) <- c("acclist", "gt1m")
	}

return(out)
}

gt3xAccDir <- function(accFileList, save, compress = "gzip", compression_level = 6){

FILES <- accFileList$FILES
FILEID <- accFileList$FILEID
path <- accFileList$path
model <- accFileList$model
N <- length(FILES)
newEnv <- new.env()
tz <- accFileList$tz

# determine saving mode

	if(is.logical(save)){
		if(save){
			saveDir <- paste(path, "/accRfiles", sep = "")
			if(file.exists(saveDir)){
				cat(paste("Directory ", saveDir, " already exists. \n Do you still want to save R files in that folder (y/n)? \n", sep = ""))
				ans <- readLines(stdin(),1)
				if(as.character(ans) == "n")
					stop("Data processing stopped by user")
						else if(as.character(ans) != "y") stop("Ambiguous answer. Process halted.")
				} else {
				dir.create(saveDir)
				cat("Directory ", saveDir, " created \n")
				}
		} else {saveDir <- NULL}
	}

	if(is.character(save)){
		saveDir <- save
		if(file.exists(saveDir)){
			cat(paste("Directory ", saveDir, " already exists. \n Do you still want to save R files in that folder (y/n)? \n", sep = ""))
			ans <- readLines(stdin(),1)
				if(as.character(ans) == "n")
					stop("Data processing stopped by user")
				else if(as.character(ans) != "y") stop("Ambiguous answer. Process halted.")
		}
		if(!file.exists(saveDir)){
			cat("Directory for output does not exist. \n Do you want to create it (y/n)? \n")
			ans <- readLines(stdin(),1)
				if(as.character(ans) == "n")
					stop("Data processing stopped by user")
				else if(as.character(ans) != "y") stop("Ambiguous answer. Process halted.")
					else if(as.character(ans) == "y"){
						dir.create(saveDir)
						cat("Directory ", saveDir, " created \n")
					}
		}
	}

# Loop over files

summaryFile <- data.frame(fileid = "", serial = "", nobs = 0, epoch = 0, mode = NA, ts_start = as.POSIXlt("2000/01/01 00:00:00"), tz = "", voltage = NA, ts_dl = as.POSIXlt("2000/01/01 00:00:00"), stringsAsFactors = FALSE)

pb <- winProgressBar(title = "progress bar", min = 0, max = N, width = 300)
	if(!is.null(saveDir)){
		for(i in 1:N){
			object <- gt3xAccFile(file = FILES[i], path = path, fileid = FILEID[i], tz = tz, sparse = accFileList$sparse)
			summaryFile[i,] <- object$info
			assign(x = FILEID[i], value = object, envir = newEnv)
			save(list = c(FILEID[i]), file = paste(saveDir, "/", FILEID[i], ".Rdata", sep =""), envir = newEnv, compress = compress, compression_level = compression_level)
			remove(list = c(FILEID[i]), envir = newEnv)
			setWinProgressBar(pb, i, title = paste(round(i/N*100, 0), "% done"))
		}
		close(pb)
		cat("Output is ready \n")
		out <- summaryFile
		class(out) <- c("acclist", "gt3x", "summaryFile")
	} else {
		out <- vector("list", N)
		for(i in 1:N){
			out[[i]] <- gt3xAccFile(file = FILES[i], path = path, fileid = FILEID[i], tz = tz, sparse = accFileList$sparse)
			summaryFile[i,] <- out[[i]]$info
			setWinProgressBar(pb, i, title = paste(round(i/N*100, 0), "% done"))
		}
		close(pb)
		names(out) <- FILEID
		cat("Output is ready \n")
		attr(out, "info") <- summaryFile
		class(out) <- c("acclist", "gt3x")
	}

return(out)

}

gt1mAccFile <- function(file, path, fileid, counts.pos = 1, tz = "Europe/London", sparse = FALSE, fault = 32767){

filename <- paste(path, file, sep = "/")

fileConnection <- file(filename, "r")
	if(isOpen(fileConnection, "r")) {
		seek(fileConnection, 0, rw = "r")
		Lines <- readLines(fileConnection)
		close(fileConnection)
	}

# Determine device mode	and epoch
Mode <- grep("Mode", Lines)
Mode <- as.numeric(strsplit(Lines[Mode], "= ")[[1]][2])
epoch <- grep("Epoch", Lines)
epoch <- strsplit(Lines[epoch], ") ")[[1]][2]
epoch <- as.numeric(strsplit(epoch, ":")[[1]])
epoch <- epoch[1]*3600 + epoch[2]*60 + epoch[3]
serial <- grep("Serial", Lines)
serial <- strsplit(Lines[serial], ":")[[1]][2]
serial <- sub("^\\s+", "", serial)

# Other info
Voltage <- grep("Voltage", Lines)
Voltage <- as.numeric(strsplit(strsplit(Lines[Voltage], "Mode")[[1]][1], ":")[[1]][2])

# Download date and time
sel <- grep("Download Time", Lines)
downTime <- gsub("Download Time ", "", Lines[sel])
downTime <- gsub("[[:blank:]]", "", downTime)
sel <- grep("Download Date ", Lines)
downDate <- gsub("Download Date ", "", Lines[sel])
downDate <- gsub("[[:blank:]]", "", downDate)
downDate <- strsplit(downDate, "/")[[1]]

	if (nchar(downDate[1]) == 1) {
		downDate[1] <- paste("0", downDate[1], sep = "")
	}
	if (nchar(downDate[2]) == 1) {
		downDate[2] <- paste("0", downDate[2], sep = "")
	}

downDate <- infoDate(downDate, fileid)$date
TS_dl <- paste(downDate, downTime, sep = " ")
TS_dl <- as.POSIXlt(TS_dl, tz = tz)

# Start date and time
sel <- grep("Start Time", Lines)
startTime <- gsub("Start Time ", "", Lines[sel])
startTime <- gsub("[[:blank:]]", "", startTime)
sel <- grep("Start Date ", Lines)
startDate <- gsub("Start Date ", "", Lines[sel])
startDate <- gsub("[[:blank:]]", "", startDate)
startDate <- strsplit(startDate, "/")[[1]]

	if (nchar(startDate[1]) == 1) {
		startDate[1] <- paste("0", startDate[1], sep = "")
	}
	if (nchar(startDate[2]) == 1) {
		startDate[2] <- paste("0", startDate[2], sep = "")
	}
	
startDate <- infoDate(startDate, fileid)	
TS_orig <- paste(startDate$date, startTime, sep = " ")
TS_orig <- as.POSIXlt(TS_orig, tz = tz)
tz <- format(TS_orig, "%Z")
	if(!(tz %in% c("GMT","BST")))
		{warning(paste("GMT/BST not determined for accelerometer ", fileid)); tz <- "NA"}

# Process counts and steps

startL <- grep("-----", Lines)[2] + 1
endL <- length(Lines)
tmp <- Lines[startL:endL]
tmp <- lapply(tmp, function(x) strsplit(gsub("[[:blank:]]+", " ", x), " ")[[1]])
accData <- as.numeric(unlist(lapply(tmp, function(x) x[x!=""])))
N <- length(accData)

# Check consistency data

oddN <- FALSE
	if(Mode == 1){
		if(!(counts.pos %in% (1:2)))
			stop("'counts.pos' is either 1 or 2")
		steps.pos <- if(counts.pos == 1) 2 else 1
		cts <- accData[seq(counts.pos, N, by = 2)]
		steps <- accData[seq(steps.pos, N, by = 2)]
		n <- length(cts)
		if(N%%2 != 0){
			paste("Numbers of counts and steps do not match for accelerometer ", fileid,
			". Adding '0' as last value for shorter vector", sep = "");
			oddN <- TRUE;
			n1 <- length(cts);
			n2 <- length(steps);
			if(n1 < n2) {cts <- c(cts,0)} else {steps <- c(steps,0)};
			n <- max(c(n1,n2))
			}
		error_c <- errorChk(cts, fault = fault)
		error_s <- errorChk(steps, fault = fault)
	} else if(Mode == 0){
		cts <- accData
		steps <- 0
		n <- N
		error_c <- errorChk(cts, fault = fault)
		error_s <- 0
	}

# Summary file
	
summaryFile <- data.frame(fileid = fileid, serial = serial, nobs = n, epoch = epoch, mode = Mode, ts_start = TS_orig, tz = tz, voltage = Voltage, ts_dl = TS_dl, stringsAsFactors = FALSE)

# Store in sparse matrix format
err_s <- list(fileid = fileid, counts = table(error_c), steps = table(error_s), date = startDate$date_code, odd_number = oddN)

Data <- data.frame(counts = cts, steps = steps, error_c = error_c, error_s = error_s)

out <- if(sparse){
	list(df = as.matrix.csr(as.matrix(Data)), info = summaryFile, error_summary = err_s)
	} else {
	list(df = Data, info = summaryFile, error_summary = err_s)
	}

attr(out, "sparse") <- sparse
attr(out, "labels") <- colnames(Data)
class(out) <- c("accfile", "gt1m")
return(out)

}

gt3xAccFile <- function(file, path, fileid, tz = "Europe/London", sparse = FALSE, fault = 32767){

filename <- paste(path, file, sep = "/")

fileConnection <- file(filename, "r")
	if(isOpen(fileConnection, "r")) {
		seek(fileConnection, 0, rw = "r")
		Lines <- readLines(fileConnection)
		close(fileConnection)
	}

# Determine device mode	and epoch
Mode <- grep("Mode", Lines)
Mode <- as.numeric(strsplit(Lines[Mode], "= ")[[1]][2])
epoch <- grep("Epoch", Lines)
epoch <- strsplit(Lines[epoch], ") ")[[1]][2]
epoch <- as.numeric(strsplit(epoch, ":")[[1]])
epoch <- epoch[1]*3600 + epoch[2]*60 + epoch[3]
serial <- grep("Serial", Lines)
serial <- strsplit(Lines[serial], ":")[[1]][2]
serial <- sub("^\\s+", "", serial)

# Other info
Voltage <- grep("Voltage", Lines)
Voltage <- as.numeric(strsplit(strsplit(Lines[Voltage], "Mode")[[1]][1], ":")[[1]][2])
formatDate <- strsplit(strsplit(Lines[1], "date format")[[1]][2], "at")[[1]][1]

if(is.na(formatDate)){


# Download date and time
sel <- grep("Download Time", Lines)
downTime <- gsub("Download Time ", "", Lines[sel])
downTime <- gsub("[[:blank:]]", "", downTime)
sel <- grep("Download Date ", Lines)
downDate <- gsub("Download Date ", "", Lines[sel])
downDate <- gsub("[[:blank:]]", "", downDate)
downDate <- strsplit(downDate, "/")[[1]]

	if (nchar(downDate[1]) == 1) {
		downDate[1] <- paste("0", downDate[1], sep = "")
	}
	if (nchar(downDate[2]) == 1) {
		downDate[2] <- paste("0", downDate[2], sep = "")
	}

downDate <- infoDate(downDate, fileid)$date
TS_dl <- paste(downDate, downTime, sep = " ")
TS_dl <- as.POSIXlt(TS_dl, tz = tz)

# Start date and time
sel <- grep("Start Time", Lines)
startTime <- gsub("Start Time ", "", Lines[sel])
startTime <- gsub("[[:blank:]]", "", startTime)
sel <- grep("Start Date ", Lines)
startDate <- gsub("Start Date ", "", Lines[sel])
startDate <- gsub("[[:blank:]]", "", startDate)
startDate <- strsplit(startDate, "/")[[1]]

	if (nchar(startDate[1]) == 1) {
		startDate[1] <- paste("0", startDate[1], sep = "")
	}
	if (nchar(startDate[2]) == 1) {
		startDate[2] <- paste("0", startDate[2], sep = "")
	}
	
startDate <- infoDate(startDate, fileid)	
TS_orig <- paste(startDate$date, startTime, sep = " ")
TS_orig <- as.POSIXlt(TS_orig, tz = tz)

} else {

formatDate <- gsub("^\\s+|\\s+$", "", formatDate)
formatDate <- strsplit(formatDate, "/")[[1]]
formatDate[formatDate=="M"] <- "m"
formatDate[formatDate=="D"] <- "d"
formatDate[formatDate=="yyyy"] <- "Y"
formatDate <- paste("%", formatDate, collapse = "/", sep = "")

# Download date and time
sel <- grep("Download Time", Lines)
downTime <- gsub("Download Time ", "", Lines[sel])
downTime <- gsub("[[:blank:]]", "", downTime)
sel <- grep("Download Date ", Lines)
downDate <- gsub("Download Date ", "", Lines[sel])
downDate <- gsub("[[:blank:]]", "", downDate)

TS_dl <- strptime(paste(downDate,downTime), format = paste(formatDate, " %H:%M:%S", sep = ""), tz = tz)

# Determine start date and time

sel <- grep("Start Time", Lines)
startTime <- gsub("Start Time ", "", Lines[sel])
startTime <- gsub("[[:blank:]]", "", startTime)
sel <- grep("Start Date ", Lines)
startDate <- gsub("Start Date ", "", Lines[sel])
startDate <- gsub("[[:blank:]]", "", startDate)
startDate <- strsplit(startDate, "/")[[1]]

	if (nchar(startDate[1]) == 1) {
		startDate[1] <- paste("0", startDate[1], sep = "")
	}
	if (nchar(startDate[2]) == 1) {
		startDate[2] <- paste("0", startDate[2], sep = "")
	}

startDate <- infoDate(x = startDate, id = fileid, format = formatDate)

TS_orig <- strptime(paste(startDate$date,startTime), format = "%Y-%d-%m %H:%M:%S", tz = tz)
}

tz <- format(TS_orig, "%Z")
	if(!(tz %in% c("GMT","BST")))
		{warning(paste("GMT/BST not determined for accelerometer ", fileid)); tz <- "NA"}

# Process data

startL <- grep("-----", Lines)[2] + 1
endL <- length(Lines)
tmp <- Lines[startL:endL]
tmp <- lapply(tmp, function(x) strsplit(gsub("[[:blank:]]+", " ", x), " ")[[1]])
ncols <- length(strsplit(tmp[[1]], ",")[[1]])
if(ncols > 4) warning(paste("Number of accelerometer variables is", ncols))
accData <- matrix(as.numeric(unlist(strsplit(unlist(tmp), ","))), ncol = ncols, byrow = TRUE)
n <- nrow(accData)
colnames(accData) <- c('y','x','z','steps')[1:ncols]
error <- NULL
for(j in 1:ncols){
error <- cbind(error, errorChk(accData[,j], fault = fault))
}
colnames(error) <- paste("error", substr(colnames(accData), 1, 1), sep = "_")

# Summary file
	
summaryFile <- data.frame(fileid = fileid, serial = serial, nobs = n, epoch = epoch, mode = Mode, ts_start = TS_orig, tz = tz, voltage = Voltage, ts_dl = TS_dl, stringsAsFactors = FALSE)

# Error summary
err_s <- apply(error, 2, table)
err_s <- c(fileid = fileid, as.list(err_s), date = startDate$date_code, odd_number = NA)

# Output

Data <- data.frame(accData, error)

out <- if(sparse){
	list(df = as.matrix.csr(Data), info = summaryFile, error_summary = err_s)
} else {
	list(df = data.frame(Data), info = summaryFile, error_summary = err_s)
	}

attr(out, "sparse") <- sparse
attr(out, "labels") <- colnames(accData)
class(out) <- c("accfile", "gt3x")
return(out)


}

infoDate <- function(x, id, ...){

# 0 y/m/d
# 1 y/d/m
# 2 not recognised

y <- paste(x[3], x[2], x[1], sep = "-")
test <- try(as.POSIXlt(y, ...), silent = TRUE)
date_code <- 0
	if(class(test)[1]=="try-error"){
		y <- paste(x[3], x[1], x[2], sep = "-")
		warning(paste("Date format for accelerometer", id, "changed to '%Y-%m-%d'"))
		date_code <- date_code + 1
		test <- try(as.POSIXlt(y, ...), silent = TRUE)
		if(class(test)[1]=="try-error")
		{warning(paste("Date format for accelerometer", id, "set to '01/01/2099'")); y <- "01/01/2099";
			date_code <- date_code + 1}
	}

return(list(date = y, date_code = date_code))
}

# Print 'acclist' object

print.acclist <- function(x, ...){

y <- if("summaryFile" %in% class(x)) x else attributes(x)$info
class(y) <- "data.frame"

rd <- format(range(y$ts_start, na.rm = TRUE), "%Y-%m-%d")
n <- length(y$fileid)
	if(n > 1){
		cat(paste("There are ", n, " accelerometer files. \n", "Date range is ", rd[1], " to ", rd[2], "\n", sep = ""))
		cat("\n")
		print(y, quote = FALSE)
	} else {print(y, quote = FALSE)}
}

# Print accfile

print.accfile <- function(x, ...){

print(x$info, quote = FALSE)

}

### Classify errors

errorChk <- function(x, fault = 32767){

n <- length(x)
Error <- rep(0,n)

# NAs
NAS <- is.na(x)
if(sum(NAS) > 0){
Error[NAS] <- 3
}

# aberrant
if(all(x[!NAS] > 10000)) Error[!NAS] <- 1
if(all(x[!NAS] == min(x[!NAS]))) Error[!NAS] <- 1
if(min(x[!NAS]) > 0) Error[!NAS] <- 1

Error[x == fault & !NAS] <- 1
Error[x == fault & !NAS] <- 1

# negative
Error[x < 0 & !NAS] <- 2


return(Error)
}


# Classify wear/non-wear time

markwear.accfile <- function(object, value, which = "counts", rescale.epoch = 60, nz = 0, keep.error = FALSE){

# consecutive zero-counts (value is expressed in minutes)
sparse <- attr(object, "sparse")
info <- object$info
if(info$epoch < 1) {
	stop("Epochs less than 1 second are not allowed")
	} else {f <- rescale.epoch/info$epoch}

if(sparse){
	Data <- as.data.frame(as.matrix(object$df))
	colnames(Data) <- attr(object, "labels")
} else {
	Data <- object$df
}

nn <- intersect(c("x","y","z","counts","steps"), colnames(Data))

if("gt1m" %in% class(object)){
	if(!which %in% nn) stop(cat("Argument 'which' must be one of", nn,"\n"))
	x <- Data[,which]
	err <- paste("error", substr(which, 1, 1), sep = "_")
	err <- Data[,err]
} else if("gt3x" %in% class(object)){
	if(!which %in% nn) stop(cat("Argument 'which' must be one of", nn,"\n"))
	x <- Data[,which]
	err <- paste("error", substr(which, 1, 1), sep = "_")
	err <- Data[,err]
}

value <- value*f # 'value' is converted into seconds
nz <- nz*f # 'nz' is converted into seconds
z <- rle(x)

if(nz > 0){

# detect sequences of zeroes with lenght >= value allowing for non-zero sequences not longer than 'nz'

n <- length(z$values)
from <- 1:(n-2)
middle <- from + 1
to <- from + 2

sel.v <- z$values[from] == 0 & z$values[to] == 0 & z$values[middle] != 0
sel.l <- (z$lengths[from] + z$lengths[to]) >= value & z$lengths[middle] <= nz
sel <- sort(c(from[sel.v & sel.l], middle[sel.v & sel.l], to[sel.v & sel.l]))
sel <- sel[!duplicated(sel)]

nsel <- (1:n)[!((1:n) %in% sel)]
z$values[sel] <- 0
z$values[nsel] <- ifelse(z$values[nsel]==0 & z$lengths[nsel] >= value, 0, 1) # detect sequences of zeroes with lenght >= value and replace with zeros
z <- inverse.rle(z)
} else {
z$values <- ifelse(z$values==0 & z$lengths >= value, 0, 1) # detect sequences of zeroes with lenght >= value and replace with zeros
z <- inverse.rle(z)
}

# handle errors
z <- handleError(z, err, code = "all", na = TRUE, keep.error = keep.error)

z <- factor(z, levels = c(0,1), labels = c("Non-wear","Wear"))
return(z)

}

markwear.acclist <- function(object, value, which = "counts",  rescale.epoch = 60, nz = 0, keep.error = FALSE){

# consecutive zero-counts (value is expressed in minutes)
fileids <- attributes(object)$info$fileid
N <- length(fileids)
out <- vector("list", N)

pb <- winProgressBar(title = "progress bar", min = 0, max = N, width = 300)
for(i in 1:N){
	out[[i]] <- do.call(markwear.accfile, args = list(object = object[[i]], value = value, which = which, rescale.epoch = rescale.epoch, nz = nz, keep.error = keep.error))
	setWinProgressBar(pb, i, title = paste(round(i/N*100, 0), "% done"))
}
close(pb)
names(out) <- fileids
return(out)

}

# Classify physical activity 

markpa.accfile <- function(object, value, which = "counts", rescale.epoch = 60, labels = NULL, extreme = NULL, keep.error = FALSE){

sparse <- attr(object, "sparse")
info <- object$info
if(info$epoch < 1) {
	stop("Epochs less than 1 second are not allowed")
	} else {f <- rescale.epoch/info$epoch}

if(sparse) {
	Data <- as.data.frame(as.matrix(object$df))
	colnames(Data) <- attr(object, "labels")
} else {
	Data <- object$df
}

nn <- intersect(c("x", "y", "z", "counts", "steps"), colnames(Data))
if ("gt1m" %in% class(object)) {
	if (!which %in% nn) 
		stop(cat("Argument 'which' must be one of", nn, "\n"))
	x <- Data[, which]
	err <- paste("error", substr(which, 1, 1), sep = "_")
	err <- Data[, err]
}
else if ("gt3x" %in% class(object)) {
	if (!which %in% nn) 
		stop(cat("Argument 'which' must be one of", nn, "\n"))
	x <- Data[, which]
	err <- paste("error", substr(which, 1, 1), sep = "_")
	err <- Data[, err]
}

vl <- length(value)
if(value[1] != 0) stop("First breakpoint must be 0")

if(!is.null(labels)){
	nl <- length(labels)
	if(nl != vl) stop(paste("The number of labels (", nl, ") does not match the number of possible intervals from argument value (", vl, ")", sep = ""))
} else {
	labels <- paste("PA", 1:vl, sep = "")
}

# handle  errors
x <- handleError(x, err, code = "all", na = TRUE, keep.error = keep.error)

# classify values by epoch based on 'value' breaks. NB, the time interval for value (expressed in counts per x seconds) may be different from epoch (e.g., 15 seconds). Therefore, value is divided by f = rescale.epoch/epoch 

z <- findInterval(x, vec = value/f, all.inside = F)
z <- factor(z, levels = 1:vl, labels = labels)

attr(z, "extreme") <- if(is.null(extreme)) NULL else if(extreme == "last") labels[vl] else labels[extreme]

return(z)

}

markpa.acclist <- function(object, value, which = "counts", rescale.epoch = 60, labels = NULL, extreme = NULL, keep.error = FALSE){

fileids <- attributes(object)$info$fileid
N <- length(fileids)
out <- vector("list", N)

pb <- winProgressBar(title = "progress bar", min = 0, max = N, width = 300)
for(i in 1:N){
	out[[i]] <- do.call(markpa.accfile, args = list(object = object[[i]], value = value, which = which, rescale.epoch = rescale.epoch, labels = labels, extreme = extreme, keep.error = keep.error))
	setWinProgressBar(pb, i, title = paste(round(i/N*100, 0), "% done"))
}
close(pb)
names(out) <- fileids
return(out)

}

markbouts.accfile <- function(object, value, which = "counts", bts = c(0,10,20,Inf), rescale.epoch = 60, collapse.by = "%Y-%m-%d", value.labels = NULL, bouts.labels = NULL, extreme = NULL, keep.error = FALSE){

sparse <- attr(object, "sparse")
info <- object$info

if (info$epoch < 1) {
	stop("Epochs less than 1 second are not allowed")
}
else {
	f <- rescale.epoch/info$epoch
}
if (sparse) {
	Data <- as.data.frame(as.matrix(object$df))
	colnames(Data) <- attr(object, "labels")
}
else {
	Data <- object$df
}
nn <- intersect(c("x", "y", "z", "counts", "steps"), colnames(Data))
if ("gt1m" %in% class(object)) {
	if (!which %in% nn) 
		stop(cat("Argument 'which' must be one of", nn, "\n"))
	x <- Data[, which]
	err <- paste("error", substr(which, 1, 1), sep = "_")
	err <- Data[, err]
}
else if ("gt3x" %in% class(object)) {
	if (!which %in% nn) 
		stop(cat("Argument 'which' must be one of", nn, "\n"))
	x <- Data[, which]
	err <- paste("error", substr(which, 1, 1), sep = "_")
	err <- Data[, err]
}

vl <- length(value)
if(value[1] != 0) stop("First breakpoint must be 0")

if(!is.null(value.labels)){
nl <- length(value.labels)
	if(nl != vl) stop(paste("The number of labels (", nl, ") does not match the number of possible intervals from argument value (", vl, ")", sep = ""))
} else {
value.labels <- paste("PA", 1:vl, sep = "")
}


# collapse factor
timestamp <- tsFormat(object)
collapse.factor <- as.character(format(timestamp, collapse.by))
ucf <- unique(collapse.factor)
nucf <- length(ucf)

# handle  errors
x <- handleError(x, err, code = "all", na = TRUE, keep.error = keep.error)

# classify values by epoch based on 'value' breaks. NB, the time interval for value (expressed in counts per x seconds) may be different from epoch (e.g., 15 seconds). Therefore, value is divided by f = rescale.epoch/epoch 

z <- findInterval(x, vec = value/f, all.inside = F)
uz <- sort(unique(z))
bouts <- rle(z)
bouts$lengths <- bouts$lengths/f
bouts <- xtabs(~bouts$lengths+bouts$values) # frequency table bouts by duration and PA intensity
duration <- as.numeric(rownames(bouts)) # duration of bouts
duration_factor <- cut(duration, breaks = bts, labels = bouts.labels, include.lowest = TRUE, right = FALSE)

res <- array(0, dim = c(length(bts) - 1, vl, 3, nucf), dimnames = list(levels(duration_factor), 1:vl, c("tot_duration","frequency","mean_duration"), ucf))

for(i in 1:nucf){
	x.sub <- x[collapse.factor == ucf[i]]
	z <- findInterval(x.sub, vec = value/f, all.inside = F)
	bouts <- rle(z)
	bouts$lengths <- bouts$lengths/f

	bouts <- xtabs(~bouts$lengths+bouts$values) # frequency table bouts by duration and PA intensity
	duration <- as.numeric(rownames(bouts)) # duration of bouts
	duration_factor <- cut(duration, breaks = bts, labels = bouts.labels, include.lowest = TRUE, right = FALSE)
	duration_bouts <- sweep(bouts, 1, duration, "*") # time spent in bouts
	group_bouts_dur <- apply(duration_bouts, 2, function(x,y) tapply(x, list(y), sum, na.rm = TRUE), y = duration_factor)
	if(sum(is.na(group_bouts_dur)) > 0){
		group_bouts_dur[is.na(group_bouts_dur)] <- 0
		warning("Some cells are empty. Bouts duration set to 0")
	}
	group_bouts_freq <- apply(bouts, 2, function(x,y) tapply(x, list(y), sum, na.rm = TRUE), y = duration_factor)
	if(sum(is.na(group_bouts_freq)) > 0){
		group_bouts_freq[is.na(group_bouts_freq)] <- 0
		warning("Some cells are empty. Bouts frequence set to 0")
	}
	group_bouts_meandur <- group_bouts_dur/group_bouts_freq
	
	if(!is.null(dim(group_bouts_dur))){
		res[,match(colnames(group_bouts_dur), uz),1,i] <- group_bouts_dur
		res[,match(colnames(group_bouts_freq), uz),2,i] <- group_bouts_freq
		res[,match(colnames(group_bouts_meandur), uz),3,i] <- group_bouts_meandur
	} else {
		res[,match(names(group_bouts_dur), uz),1,i] <- group_bouts_dur
		res[,match(names(group_bouts_freq), uz),2,i] <- group_bouts_freq
		res[,match(names(group_bouts_meandur), uz),3,i] <- group_bouts_meandur	
	}
}

dimnames(res)[[2]] <- value.labels
attr(res, "extreme") <- if(is.null(extreme)) NULL else if(extreme == "last") value.labels[vl] else value.labels[extreme]

return(res)

}

markbouts.acclist <- function(object, value, which = "counts", bts = c(0,10,20,Inf), rescale.epoch = 60, collapse.by = "%Y-%m-%d", value.labels = NULL, bouts.labels = NULL, extreme = NULL, keep.error = FALSE){

fileids <- attributes(object)$info$fileid
N <- length(fileids)
out <- vector("list", N)

pb <- winProgressBar(title = "progress bar", min = 0, max = N, width = 300)

for(i in 1:N){
	out[[i]] <- do.call(markbouts.accfile, args = list(object = object[[i]], value = value, which = which, bts = bts, rescale.epoch = rescale.epoch, collapse.by = collapse.by, value.labels = value.labels, bouts.labels = bouts.labels, extreme = extreme, keep.error = keep.error))
	setWinProgressBar(pb, i, title = paste(round(i/N*100, 0), "% done"))
}

close(pb)
names(out) <- fileids
return(out)

}

# Auxiliary functions

tsFromEpoch <- function(object, x){
	if(any(x <= 0)) stop("x must be positive")
	if(any(x > object$info$nobs)) warning("Timestamp is outside observed time interval")
		object$info$ts_start + object$info$epoch*(as.integer(x) - 1)
}

epochFromTS <- function(object, x){
	difft <- difftime(x, object$info$ts_start, units = "secs")
	if(difft < 0) stop(paste("Specify timestamp not prior to ", object$info$ts_start, sep =""))
	val <- abs(difft/object$info$epoch) + 1
	as.numeric(floor(val))
}


tsFormat <- function(object){

info <- object$info

# total seconds
N <- info$nobs*info$epoch

# time stamping (it accounts for BST)
TimeStamp <- info$ts_start + seq(1,N,by=info$epoch)

return(TimeStamp)
}


dateSummary <- function(object, wear, timestamp, minval  = 0, rescale.epoch = 60, keep.error = FALSE){

info <- object$info

# total seconds
N <- info$nobs*info$epoch

# start and end dates
start <- info$ts_start
end <- start + N

# interval in days
int <- as.numeric(format(end, "%j")) - as.numeric(format(start, "%j")) + 1

# days and hours
days <- as.character(format(timestamp, "%Y-%m-%d"))
hrs <- as.character(format(timestamp, "%H"))

td <- data.frame(fileid = info$fileid, table(days))
names(td)[names(td) == "Freq"] <- "freq"
td$hour_day <- td$freq*info$epoch/3600

tmp <- tapply(hrs, list(day = days), function(x) c(min(x),max(x)))
td[,c("start_day","end_day")] <- t(sapply(tmp, function(x) x))

# valid day and ID (minval is in minutes)

td$valid_mins <- tapply(as.numeric(wear)-1, list(day = days), function(x,n,rescale.epoch) sum(x, na.rm = TRUE)*n/rescale.epoch, n = info$epoch, rescale.epoch = rescale.epoch)

vd <- as.numeric(td$valid_mins >= minval)
td$IsEndDate <- td$IsStartDate <- 0
a <- seq(along.with = vd)

sels <- a[cumsum(vd) == 1][1]
if(sum(vd) != 0) td$IsStartDate[sels] <- 1
sele <- a[cumsum(vd) == sum(vd)][1]
if(sum(vd) != 0) td$IsEndDate[sele] <- 1
td$IsTruncated <- if(sum(vd) != 0) as.numeric(a < sels | a > sele) else 1

return(td)
}


# Aggregate accelerometer data

aggAccFile <- function(object, by, which = "counts", x = NULL, keep.error = FALSE){

info <- object$info
sparse <- attr(object, "sparse")

if(info$epoch > by) stop(paste("Epoch is longer than ", by, sep = ""))
f <- by/info$epoch
#if(f == 1) return(object)
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
if(!is.wholenumber(f)) stop("Only multiples of epoch length are allowed")

if (sparse) {
	Data <- as.data.frame(as.matrix(object$df))
	colnames(Data) <- attr(object, "labels")
}
else {
	Data <- object$df
}

if(is.null(x)){
nn <- intersect(c("x", "y", "z", "counts", "steps"), colnames(Data))
if ("gt1m" %in% class(object)) {
	if (!which %in% nn) 
		stop(cat("Argument 'which' must be one of", nn, "\n"))
	x <- Data[, which]
	err <- paste("error", substr(which, 1, 1), sep = "_")
	err <- Data[, err]
}
else if ("gt3x" %in% class(object)) {
	if (!which %in% nn) 
		stop(cat("Argument 'which' must be one of", nn, "\n"))
	x <- Data[, which]
	err <- paste("error", substr(which, 1, 1), sep = "_")
	err <- Data[, err]
}

# handle errors
x <- handleError(x, err, code = "all", na = TRUE, keep.error = keep.error)

}

minn <- seq(1, info$nobs, by = f)
maxn <- seq(f, info$nobs, by = f)
if(info$nobs %% f !=0) maxn <- c(maxn, info$nobs)
if(length(minn) != length(maxn)) stop("Check 'minn' and 'maxn'")

fun.do <- function(a, b, x) sum(x[a:b], na.rm = TRUE)
x <- mapply(fun.do, a = minn, b = maxn, MoreArgs = list(x = x)) # aggregated

if(sparse) {
	x <- as.matrix.csr(x)
}

TimeStamp <- tsFromEpoch(object, minn)

out <- list(outcome = x, ts_agg = TimeStamp)
attr(out, "sparse") <- sparse
class(out) <- "accfile_agg"
return(out)
}

### HANDLE ERRORS

handleError <- function(x, err, code = "all", replace = NULL, na = TRUE, keep.error = FALSE){

if(keep.error) return(x)
if(!na & is.null(replace)) stop("'na = FALSE'. The argument 'replace' is required")

sel <- if(code == "all") err != 0 else err == code
if(sum(sel) > 0){
	x[sel] <- if(na) NA else replace
}
if(sum(sel) > 0 & na) warning("NAs imputed where errors found")

return(x)

}

# Functions for data collapsing

fun.collapse <- function(x, fun = list(mean = function(x) mean(x, na.rm = TRUE), median = function(x) median(x, na.rm = TRUE),
	sd = function(x) sd(x, na.rm = TRUE))){

n <- length(fun)

if(n == 1){
	out <- match.fun(fun[[1]])(x)
	names(out) <- names(fun)
} else {
	out <- rep(NA, n)
	for(i in 1:n){
			val <- match.fun(fun[[i]])(x)
			out[i] <- if(length(val) != 1) NA else val
	}
	names(out) <- names(fun)
}

return(out)
	
}

# Collapse data

collapse.accfile <- function(object, which = "counts", palist = list(value = c(0, 100, 1e3, 5e3, 13e3), rescale.epoch = 60, labels = NULL, extreme = NULL), mwlist = list(value = 20, nz = 0, rescale.epoch = 60), collapse.by = "%Y-%m-%d", collapse.epoch = 60, aggregate.by = NULL, FUN.list = list(mean = function(x) mean(x, na.rm = TRUE)), keep.extreme = FALSE, keep.error = FALSE, ...){

Call <- match.call()

sparse <- attr(object, "sparse")
info <- object$info
if (info$epoch < 1) {
	stop("Epochs less than 1 second are not allowed")
}
if (sparse) {
	Data <- as.data.frame(as.matrix(object$df))
	colnames(Data) <- attr(object, "labels")
}
else {
	Data <- object$df
}
nn <- intersect(c("x", "y", "z", "counts", "steps"), colnames(Data))
if ("gt1m" %in% class(object)) {
	if (!which %in% nn) 
		stop(cat("Argument 'which' must be one of", nn, "\n"))
	x <- Data[, which]
	err <- paste("error", substr(which, 1, 1), sep = "_")
	err <- Data[, err]
}
else if ("gt3x" %in% class(object)) {
	if (!which %in% nn) 
		stop(cat("Argument 'which' must be one of", nn, "\n"))
	x <- Data[, which]
	err <- paste("error", substr(which, 1, 1), sep = "_")
	err <- Data[, err]
}

# Classify wear-time and PA type
wear <- do.call(markwear.accfile, args = c(list(object = object, which = which, keep.error = keep.error), mwlist))
timestamp <- tsFormat(object)
DS <- dateSummary(object, wear, timestamp, keep.error = keep.error, ...)
PA_type <- do.call(markpa.accfile, args = c(list(object = object, which = which, keep.error = keep.error), palist))
upa <- unique(PA_type)

# handle errors and extreme values
x <- handleError(x, err, code = "all", na = TRUE, keep.error = keep.error)

if(!is.null(attr(PA_type, "extreme"))){
	if(any(PA_type == attr(PA_type, "extreme")) & !keep.extreme){
		x[PA_type == attr(PA_type, "extreme")] <- NA
		if(sum(PA_type == attr(PA_type, "extreme"), na.rm = TRUE) > 0) warning("NAs imputed where extreme counts found")
	}
}

# handle wear/non-wear time
x[wear == "Non-wear"] <- NA
collapse.factor <- as.character(format(timestamp, collapse.by))

# aggregate before collapsing
if(!is.null(aggregate.by)){
	object.agg <- aggAccFile(object, by = aggregate.by, x = x)
	timestamp.agg <- object.agg$ts_agg

	# check collapsing factors
	collapse.factor.agg <- as.character(format(timestamp.agg, collapse.by))
	chk <- length(setdiff(collapse.factor, collapse.factor.agg)) == 0
	if(!chk) stop("Collapsing factors differ for original and aggregated values")
	
	out.agg <- data.frame(collapse.by = collapse.factor.agg, outcome = as.matrix(object.agg$outcome), stringsAsFactors = F)
}

# Create dataframes, drop truncated days
sel <- DS$days[DS$IsTruncated == 0]
out <- data.frame(wear = wear, collapse.by = collapse.factor, PA_type = PA_type, outcome = x, stringsAsFactors = F)
levels(out$PA_type) <- c(levels(out$PA_type), "non-wear")
out$PA_type[out$wear == "Non-wear"] <- "non-wear" # create level non-wear

# Summary

pa.summary <- aggregate(out$PA_type, list(collapse.by = out$collapse.by, PA_type = out$PA_type), function(x,n,collapse.epoch) sum(!is.na(x))*n/collapse.epoch, n=info$epoch, collapse.epoch=collapse.epoch)
pa.summary <- xtabs(x ~ collapse.by + PA_type, data = pa.summary)
nnr <- dimnames(pa.summary)[[1]]
nnc <- dimnames(pa.summary)[[2]]
pa.summary <- data.frame(matrix(pa.summary, nrow = length(nnr)))
names(pa.summary) <- nnc
pa.summary$collapse.by <- nnr

if(!is.null(aggregate.by)){
	outcome.summary <- aggregate(out.agg$outcome, list(collapse.by = out.agg$collapse.by), function(x, args) fun.collapse(x, args), args = FUN.list)
} else {
	outcome.summary <- aggregate(out$outcome, list(collapse.by = out$collapse.by), function(x, args) fun.collapse(x, args), args = FUN.list)
}

# Put in matrix format and drop dates where null

nnf <- names(FUN.list)

if(length(nnf) != 1) {
colnames(outcome.summary$x) <- paste("outcome", colnames(outcome.summary$x), sep = ".")
}
# Produce final dataframe

MAT <- data.frame(fileid = as.character(info$fileid), collapse.by = outcome.summary$collapse.by)

if(length(nnf) == 1) {
MAT <- cbind(MAT, outcome.summary[,2])
names(MAT)[3] <- paste("outcome", nnf, sep = ".")
} else {
MAT <- cbind(MAT, outcome.summary$x)
}

val <- list(outcome = merge(MAT, pa.summary))
val$call <- Call
class(val) <- "accfile.collapse"
return(val)

}


# Plot functions

plot.gt1m <- function(x, y = NULL, xlab, ylab, main, keep.error = TRUE, which = "counts", select = 1,...){

select <- select[1]
if(!which %in% c("counts","steps")) stop("Argument 'which' must be one of  c('counts','steps')")

object <- x
if("accfile" %in% class(object)){

	sparse <- attr(object, "sparse")

	if(sparse){
		x <- as.numeric(as.matrix(object$counts))
		y <- as.numeric(as.matrix(object$steps))
		errs_c <- as.numeric(as.matrix(object$error_c))
		errs_s <- as.numeric(as.matrix(object$error_s))
	} else {
		x <- as.numeric(object$df$counts)
		y <- as.numeric(object$df$steps)	
		errs_c <- as.numeric(object$df$error_c)
		errs_s <- as.numeric(object$df$error_s)
	}


	# handle errors
	x <- handleError(x, errs_c, code = "all", na = TRUE, keep.error = keep.error)
	y <- handleError(y, errs_s, code = "all", na = TRUE, keep.error = keep.error)

	z <- switch(which,
		counts = x,
		steps = y)

	timestamp <- tsFormat(object)
	
	if(missing(xlab)) xlab <- "Timestamp"
	if(missing(ylab)) ylab <- which
	if(missing(main)) main <- object$info$fileid

	plot.default(timestamp, z, type = "l", xlab = xlab, ylab = ylab, main = main, ...)
} else {

	object <- object[[select]]
	sparse <- attr(object, "sparse")

	if(sparse){
		x <- as.numeric(as.matrix(object$counts))
		y <- as.numeric(as.matrix(object$steps))
		errs_c <- as.numeric(as.matrix(object$error_c))
		errs_s <- as.numeric(as.matrix(object$error_s))
	} else {
		x <- as.numeric(object$df$counts)
		y <- as.numeric(object$df$steps)	
		errs_c <- as.numeric(object$df$error_c)
		errs_s <- as.numeric(object$df$error_s)
	}


	# handle errors
	x <- handleError(x, errs_c, code = "all", na = TRUE, keep.error = keep.error)
	y <- handleError(y, errs_s, code = "all", na = TRUE, keep.error = keep.error)

	z <- switch(which,
		counts = x,
		steps = y)

	timestamp <- tsFormat(object)

	if(missing(xlab)) xlab <- "Timestamp"
	if(missing(ylab)) ylab <- which
	if(missing(main)) main <- object$info$fileid

	plot.default(timestamp, z, type = "l", xlab = xlab, ylab = ylab, main = main, ...)


}


}

plot.gt3x <- function(x, y = NULL, xlab, ylab, main, keep.error = TRUE, which = "x", select = 1,...){

select <- select[1]
if(!which %in% c("x","y","z","steps")) stop("Argument 'which' must be one of  c('x','y','z','steps')")

object <- x
if("accfile" %in% class(object)){

	sparse <- attr(object, "sparse")

	if(sparse){
		Data <- as.matrix(object$Data)
		colnames(Data) <- attr(object, "labels")
	} else {
		Data <- object$df
	}

	ncols <- ncol(Data)
	if(ncols < 5 & which == "steps") stop("Steps variable not available")
	z <- switch(which,
		x = Data$x,
		y = Data$y,
		z = Data$z,
		steps = Data$steps)

	timestamp <- tsFormat(object)
	
	if(missing(xlab)) xlab <- "Timestamp"
	if(missing(ylab)) ylab <- which
	if(missing(main)) main <- object$info$fileid

	plot.default(timestamp, z, type = "l", xlab = xlab, ylab = ylab, main = main, ...)
} else {

	object <- object[[select]]
	sparse <- attr(object, "sparse")

	if(sparse){
		Data <- as.matrix(object$df)
		colnames(Data) <- attr(object, "labels")
	} else {
		Data <- object$df
	}


	ncols <- ncol(Data)
	if(ncols < 5 & which == "steps") stop("Steps variable not available")
	z <- switch(which,
		x = Data$x,
		y = Data$y,
		z = Data$z,
		steps = Data$steps)

	timestamp <- tsFormat(object)

	if(missing(xlab)) xlab <- "Timestamp"
	if(missing(ylab)) ylab <- which
	if(missing(main)) main <- object$info$fileid

	plot.default(timestamp, z, type = "l", xlab = xlab, ylab = ylab, main = main, ...)


}


}


