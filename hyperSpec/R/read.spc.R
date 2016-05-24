### read.spc - Import Thermo Galactic's .spc file format into an hyperSpec Object
###
### C. Beleites 2009/11/29
###
#####################################################################################################

## Define constants ---------------------------------------------------------------------------------

.nul <- as.raw (0)

## header sizes 
.spc.size <- c (hdr = 512, subhdr = 32, subfiledir = 12, loghdr = 64)

.spc.default.keys.hdr2data <- c('fexper', 'fres', 'fsource')
.spc.default.keys.hdr2log  <- c('fdate', 'fpeakpt')
.spc.default.keys.log2data <- FALSE
.spc.default.keys.log2log  <- TRUE

## axis labeling ------------------------------------------------------------------------------------

## x-axis units .....................................................................................
.spc.FXTYPE <- c (expression (`/` (x, "a. u.")),                      #0
		expression (`/` (tilde (nu), cm^-1)),
		expression (`/` (lambda, (mu * m))),
		expression (`/` (lambda, nm)),
		expression (`/` (t, s)),
		expression (`/` (t, min)),
		expression (`/` (f, Hz)),
		expression (`/` (f, kHz)),
		expression (`/` (f, MHz)),
		expression (`/` (frac (m, z), frac (u, e))),        
		expression (`/` (delta, ppm)),                      # 10
		expression (`/` (t, d)),
		expression (`/` (t, a)),
		expression (`/` (Delta*tilde (nu), cm^-1)),
		expression (`/` (E, eV)),
		NA, # old version file uses label in gcatxt
		'Diode No',
		'Channel',
		expression (`/` (x, degree)),
		expression (`/` (T, degree*F)),
		expression (`/` (T, degree*C)),                     # 20
		expression (`/` (T, K)),
		'Data Point',
		expression (`/` (t, ms)),
		expression (`/` (t, micro*s)),
		expression (`/` (t, ns)),
		expression (`/` (f, GHz)),
		expression (`/` (lambda, cm)),
		expression (`/` (lambda, m)),
		expression (`/` (lambda, mm)),
		expression (`/` (t, h))                             # 30
)

.spc.xlab <- function (x) {
   if (is.character (x))
     x
   else if (x <= length (.spc.FXTYPE) + 1)
     .spc.FXTYPE [x + 1]
	else
     ## x = 255 is for double interferogram and supposed not to have a label.
     ## Thus, returning NA is appropriate
     NA
}

## y-axis units .....................................................................................
.spc.FYTYPE <- c (expression (`/` (I[Ref], "a. u.")),                      # -1
		expression (`/` (I, "a. u.")),
		expression (`/` (I[IGRM], "a. u.")),
		'A',
		expression (frac ((1 - R)^2, 2 * R)),
		'Counts',
		expression (`/` (U, V)),
		expression (`/` (y, degree)),
		expression (`/` (I, mA)),
		expression (`/` (l, mm)),
		expression (`/` (U, mV)),
		expression (-log (R)),                              # 10
		expression (`/` (y, '%')),
		expression (`/` (I, 'a. u.')),        
		expression (I / I[0]),                             
		expression (`/` (E, J)),
		NA, # old version file uses label in gcatxt
		expression (`/` (G, dB)),
		NA, # old version file uses label in gcatxt
		NA, # old version file uses label in gcatxt
		expression (`/` (T, degree*F)),
		expression (`/` (T, degree*C)),                     # 20
		expression (`/` (T, K)),
		'n',
		'K', # extinction coeaffictient
		expression (Re (y)),
		expression (Im (y)),
		'y (complex)', # complex
		'T',
		'R',
		expression (`/` (I, 'a. u.')),
		expression (`/` (I[Emission], 'a. u.'))
)
.spc.ylab <- function(x){
   if (is.character (x))
     x
   else if (x <= 26)
		.spc.FYTYPE [x + 2]
	else if (x %in% 128 : 131)
		.spc.FYTYPE [x - 99]
	else
		NA
}

## helper functions ---------------------------------------------------------------------------------
### raw.split.nul - rawToChar conversion, splitting at \0
raw.split.nul <- function (raw, trunc = c (TRUE, TRUE), firstonly = FALSE, paste.collapse = NULL) {
	# todo make better truncation
	trunc <- rep (trunc, length.out = 2)
	
	if (trunc [1] && raw [1] == .nul)
		raw <- raw [-1]
	if (trunc [2]) {
		tmp <- which (raw > .nul)
		if (length (tmp) == 0) 
			return ("")
		raw <- raw [1 : tail (tmp, 1)]	
	} 
	if (raw [length (raw)] != .nul)
		raw <- c (raw , .nul)
	
	tmp <- c (0, which (raw == .nul))
	
	out <- character (length (tmp) - 1)
	for (i in 1 : (length (tmp) - 1))
		if (tmp [i] + 1 < tmp [i + 1] - 1)
			out [i] <- rawToChar (raw [(tmp [i] + 1)  : (tmp [i + 1] - 1)])
	
	if (length (out) > 1L){
     if (firstonly){
       
       message ("multiple strings encountered in spc file ", paste (out, collapse = ", "), ": using only the first one.")
       out <- out [1]
     
       } else if (! is.null (paste.collapse)){
       
       if (hy.getOption ("debuglevel") > 2L)
         message ("multiple strings encountered in spc file ", paste (out, collapse = ", "), " => pasting.")
       
       out <- paste (out, collapse = paste.collapse)
     }
	}
	
	out
}

## file part reading functions ----------------------------------------------------------------------

## read file header .................................................................................
##
##

.spc.filehdr <- function (raw.data) {
	## check file format
	## NEW.LSB = 75 supported,
	## NEW.MSB = 76 not supported (neither by many Grams software according to spc doc)
	## OLD     = 77 not supported (replaced by new format in 1996)
	if (raw.data [2] != 75)  
		stop ("Wrong spc file format version (or no spc file at all).\n",
				"Only 'new' spc files (1996 file format) with LSB word order are supported.") 
	
	hdr <- list (ftflgs   = readBin          (raw.data [        1], "integer", 1, 1, signed = FALSE),
			## byte 2 is already interpreted
			fexper   = readBin       (raw.data [        3], "integer", 1, 1, signed = TRUE ),
			fexp     = readBin       (raw.data [        4], "integer", 1, 1, signed = TRUE ),
			fnpts    = readBin       (raw.data [  5 :   8], "integer", 1, 4                ),
			ffirst   = readBin       (raw.data [  9 :  16], "double",  1, 8                ),
			flast    = readBin       (raw.data [ 17 :  24], "double",  1, 8                ),
			fnsub    = readBin       (raw.data [ 25 :  28], "integer", 1, 4                ),
			fxtype   = readBin       (raw.data [       29], "integer", 1, 1, signed = FALSE),
			fytype   = readBin       (raw.data [       30], "integer", 1, 1, signed = FALSE),
			fztype   = readBin       (raw.data [       31], "integer", 1, 1, signed = FALSE),
			fpost    = readBin       (raw.data [       32], "integer", 1, 1, signed = TRUE ),
			fdate    = readBin       (raw.data [ 33 :  36], "integer", 1, 4                ),
			fres     = raw.split.nul (raw.data [ 37 :  45], paste.collapse = "\r\n"),
			fsource  = raw.split.nul (raw.data [ 46 :  54], paste.collapse = "\r\n"),
			fpeakpt  = readBin       (raw.data [ 55 :  56], "integer", 1, 2, signed = FALSE),
			fspare   = readBin       (raw.data [ 57 :  88], "numeric", 8, 4                ),
			fcmnt    = raw.split.nul (raw.data [ 89 : 218], paste.collapse = "\r\n"),
			fcatxt   = raw.split.nul (raw.data [219 : 248], trunc = c (FALSE, TRUE)        ),                
			flogoff  = readBin       (raw.data [249 : 252], "integer", 1, 4), #, signed = FALSE),
			fmods    = readBin       (raw.data [253 : 256], "integer", 1, 4), #, signed = FALSE),
			fprocs   = readBin       (raw.data [      257], "integer", 1, 1, signed = TRUE ),
			flevel   = readBin       (raw.data [      258], "integer", 1, 1, signed = TRUE ),
			fsampin  = readBin       (raw.data [259 : 260], "integer", 1, 2, signed = FALSE),
			ffactor  = readBin       (raw.data [261 : 264], "numeric", 1, 4                ),
			fmethod  = raw.split.nul (raw.data [265 : 312]),
			fzinc    = readBin       (raw.data [313 : 316], "numeric", 1, 4), #, signed = FALSE),
			fwplanes = readBin       (raw.data [317 : 320], "integer", 1, 4), #, signed = FALSE),
			fwinc    = readBin       (raw.data [321 : 324], "numeric", 1, 4                ),
			fwtype   = readBin       (raw.data [      325], "integer", 1, 1, signed = TRUE ),
			## 187 bytes reserved
			.last.read = .spc.size ['hdr'] 
	)
	
	## R doesn't have unsigned long int .................................
   if (any (unlist (hdr [c ("flogoff", "fmods", "fwplanes")]) < 0))
     stop ("error reading header: R does not support unsigned long integers.",
           "Please contact the maintainer of the package.")
   


   ## do some post processing ..........................................
   
	experiments <- c ("General", "Gas Chromatogram", "General Chromatogram", "HPLC Chromatogram",
			"NIR Spectrum", "UV-VIS Spectrum", "* reserved *", "X-ray diffraction spectrum",
			"Mass Spectrum", "NMR Spectrum", "Raman Spectrum", "Fluorescence Spectrum",
			"Atomic Spectrum", "Chroatography Diode Array Data")
	hdr$fexper <- factor (hdr$fexper + 1, levels = seq_along (experiments))
	levels (hdr$fexper) <- experiments
	
	hdr$ftflgs <- .spc.ftflags (hdr$ftflgs)
	
	hdr$fdate  <- ISOdate (year  = hdr$fdate %/% 1048560,
			month = hdr$fdate %/% 65536 %%  16,
			day   = hdr$fdate %/% 2048 %% 32,
			hour  = hdr$fdate %/% 64 %% 32,
			min   = hdr$fdate %% 64)
	
	## interferogram ?
	## if not, hdr$fpeakpt is set to NULL
	if (hdr$fytype == 1){ 
		if (hdr$fpeakpt != 0)
			hdr$fpeakpt <- hdr$fpeakpt + 1
	} else {
		hdr$fpeakpt <- NULL
	}
	
	## set the axis labels
	if (hdr$ftflgs ['TALABS']) {
		# TODO: find test data
		tmp <- rep (0, 4)
		tmp [seq_along (hdr$fcatxt)] <- nchar (hdr$fcatxt)
		
		if (tmp [1] > 0) hdr$fxtype <- hdr$fcatxt[1]
		if (tmp [2] > 0) hdr$fytype <- hdr$fcatxt[2]
		if (tmp [3] > 0) hdr$fztype <- hdr$fcatxt[3]
		if (tmp [4] > 0) hdr$fwtype <- hdr$fcatxt[4]
	} 
   hdr$fxtype <- .spc.xlab (hdr$fxtype)
   hdr$fytype <- .spc.ylab (hdr$fytype)
   hdr$fztype <- .spc.xlab (hdr$fztype)
   hdr$fwtype <- .spc.xlab (hdr$fwtype)
 
   
	## File with subfiles with individual x axes? 
	## Then there should be a subfile directory:
	if (hdr$ftflgs ['TXYXYS'] && hdr$ftflgs ['TMULTI']){ 
		## try to reject impossible values for the subfiledir offset
		if (hdr$fnpts > length (raw.data) || 
				(hdr$fnpts > hdr$flogoff && hdr$flogoff > 0) ||
				hdr$fnpts < 512)
			.spc.error (".spc.read.hdr", list (hdr = hdr),
					"file header flags specify TXYXYS and TMULTI, ",
					"but fnpts does not give a valid offset for the subfile directory.\n hdr$ftflgs = ",
					paste (names (hdr$ftflgs)[hdr$ftflgs], collapse = " | "),
					" (", sum (2^(0:7) [hdr$ftflgs]) , ")\n",
					"You can try to read the file using hdr$ftflgs & ! TXYXYS (",
					sum (2^(0 : 7) [hdr$ftflgs & c (TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE)]),
					"). This assumes that all subfiles do have the same x axis.\n\n")
		
		hdr$subfiledir <- hdr$fnpts
		hdr$fnpts <- 0
	} else {
		hdr$subfiledir <- 0
	}
	
	
	## some checks ......................................................
	
	if (hdr$ftflgs ['TMULTI']){
		## multiple spectra in file
		if (hdr$fnsub <= 1)
      if (hy.getOption ("debuglevel") >= 2L)
        message ("spc file header specifies multiple spectra but only zero or one subfile.")
	} else {
		## single spectrum file
		if (hdr$fnsub == 0)
			hdr$fnsub <- 1
		
		if (hdr$fnsub >  1) {
			warning ("spc file header specifies single spectrum file  but ", hdr$fnsub,
					" subfiles (spectra).\nOnly first subfile will be read.")
			hdr$fnsub <- 1
		}
		
		if (hdr$ftflgs ['TRANDM']) 
			message ("spc file header: file type flag TRANDM encountered => Enforcing TMULTI.")
		
		if (hdr$ftflgs ['TORDRD'])
			message ("spc file header: file type flag TORDRD encountered => Enforcing TMULTI.")
		
		if ((hdr$ftflgs ['TRANDM'] || hdr$ftflgs ['TORDRD']) && hdr$fnsub > 1)
			hdr$ftflgs ['TMULTI'] <- TRUE
	}
	
	if (hdr$ftflgs ['TXYXYS'] && ! hdr$ftflgs ['TXVALS']) {
		warning ("spc file header: file type flag TXYXYS encountered => Enforcing TXVALS.")
		hdr$ftflgs ['TXVALS'] <- TRUE
	}
	
	if (hdr$fwplanes > 0)
	warning ("w planes found! This is not yet tested as the developer didn't have access to such files.\n",
			"Please contact the package maintainer ", maintainer ("hyperSpec"), 
			" stating whether the file was imported successfully or not.")
	
	hdr
}

## read sub file header .............................................................................
##
## needs header for consistency checks
##

.spc.subhdr <- function (raw.data, pos, hdr) {
	subhdr <- list (subflgs   =          raw.data [pos + (      1)],
			subexp    = readBin (raw.data [pos + (      2)], "integer", 1, 1, signed = TRUE),
			subindx   = readBin (raw.data [pos + ( 3 :  4)], "integer", 1, 2, signed = FALSE),
			subtime   = readBin (raw.data [pos + ( 5 :  8)], "numeric", 1, 4),
			subnext   = readBin (raw.data [pos + ( 9 : 12)], "numeric", 1, 4),
			subnois   = readBin (raw.data [pos + (13 : 16)], "numeric", 1, 4),
			subnpts   = readBin (raw.data [pos + (17 : 20)], "integer", 1, 4), #, signed = FALSE),
			subscan   = readBin (raw.data [pos + (21 : 24)], "integer", 1, 4), #, signed = FALSE),
			subwlevel = readBin (raw.data [pos + (25 : 28)], "numeric", 1, 4))
	## 4 bytes reserved

	## R doesn't have unsigned long int .................................
   if (any (unlist (subhdr [c ("subnpts", "subscan")]) < 0))
     stop ("error reading subheader: R does not support unsigned long integers.",
           "Please contact the maintainer of the package.")
	
	hdr$.last.read <- pos + .spc.size ['subhdr']
	
	## checking
	if (subhdr$subexp == -128 && hdr$fexp != -128)
		message ("subfile ", subhdr$subindx,  " specifies data type float, but file header doesn't.",
				"\n=> Data will be interpreted as float.")
	
	if (subhdr$subnpts > 0 && subhdr$subnpts != hdr$fnpts && ! hdr$ftflgs ['TXYXYS'])
		message ('subfile ', subhdr$subindx, ": number of points in file header and subfile header ",
             "inconsistent. => Going to use subheader.")
	
	if (subhdr$subnpts == 0){
		if (hdr$ftflgs ['TXYXYS'])
			message ('subfile ', subhdr$subindx, ': number of data points per spectrum not specified. ',
					'=> Using file header information (', hdr$fnpts, ').')
		subhdr$subnpts <- hdr$fnpts
	}
	
	if (! hdr$ftflgs ['TXYXYS'])
		if (hdr$fnpts != subhdr$subnpts) {
			.spc.error (".spc.read.subhdr", list (hdr = hdr, subhdr = subhdr),
					"hdr and subhdr differ in number of points per spectrum, ",
					"but TXYXYS is not specified.\n hdr$ftflgs = ",
					paste (names (hdr$ftflgs)[hdr$ftflgs], collapse = " | "),
					" (", sum (2^(0:7) [hdr$ftflgs]) , ")\n",
					"You can try to read the file using hdr$ftflgs | TMULTI | TXYXYS (",
					sum (2^(0 : 7) [hdr$ftflgs |
											c (FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE)]),
					").\n\n")
		}
	
#  str (subhdr)
	## according to .spc file documentation:
	if (! hdr$ftflgs ['TMULTI'])
		subhdr$subexp <- hdr$fexp
   else if (hdr$fexp == -128 && subhdr$subexp != -128) {
     message ("Header file specifies float data format, but subfile uses integer exponent. ",
              "=> Using file header settings.")
     subhdr$subexp <- -128
   }
	
	## the z values
	if (hdr$fzinc == 0) # should only happen for the first subfile...
		hdr$fzinc = subhdr$subnext - subhdr$subtime
	
	if (subhdr$subindx == 0)
		hdr$firstz <- subhdr$subtime
	
	if (subhdr$subtime == 0)
		subhdr$subtime = subhdr$subindx * hdr$fzinc + hdr$firstz
	
	## the w values
	if (hdr$fwplanes > 0) {
		if (subhdr$subwlevel != 0) {
			subhdr$w <- subhdr$subwlevel
			
		} else if (subhdr$subindx %% hdr$fwplanes == 1)
			subhdr$w <- hdr$subhdr$w +  hdr$fwinc
		else
			subhdr$w <- hdr$subhdr$w
	}
	
	
	hdr$subhdr <- subhdr
	
	hdr
}

## read subfile directory ...........................................................................
##

.spc.subfiledir <- function (raw.data, pos, nsub) {
	dir <- data.frame (ssfposn = rep (NA, nsub),
			ssfsize = rep (NA, nsub),
			ssftime = rep (NA, nsub))
	
	for (s in seq_len (nsub)){
		dir [s,] <- c (readBin (raw.data [pos + ( 1 :  4)], "integer", 1, 4), # , signed = FALSE),
				readBin (raw.data [pos + ( 5 :  8)], "integer", 1, 4), # , signed = FALSE),
				readBin (raw.data [pos + ( 9 : 12)], "numeric", 1, 4))
		pos <- pos + .spc.size ['subfiledir']
	}

	## R doesn't have unsigned long int .................................
   if (any (dir [, 1:2] < 0))
     stop ("error reading subfiledir: R does not support unsigned long integers.",
           "Please contact the maintainer of the package.")
   
#	dir$ssfposn <- dir$ssfposn
	dir
}

## read log block header ............................................................................
##

.spc.log <- function (raw.data, pos, log.bin, log.disk, log.txt, keys.log2data,  keys.log2log,
                      replace.nul = as.raw (255), iconv.from = "latin1", iconv.to = "utf8") {
	
	if (pos == 0) # no log block exists
		return (list (data = list (),
						log = list ()))
	
	loghdr <- list (logsizd = readBin (raw.data [pos + ( 1 :  4)], "integer", 1, 4), #  , signed = FALSE),
			logsizm = readBin (raw.data [pos + ( 5 :  8)], "integer", 1, 4), # , signed = FALSE),
			logtxto = readBin (raw.data [pos + ( 9 : 12)], "integer", 1, 4), # , signed = FALSE),
			logbins = readBin (raw.data [pos + (13 : 16)], "integer", 1, 4), # , signed = FALSE),
			logdsks = readBin (raw.data [pos + (17 : 20)], "integer", 1, 4), # , signed = FALSE),
			## 44 bytes reserved
			.last.read = pos + .spc.size ['loghdr']
	)
	
	## R doesn't have unsigned long int .................................
   if (any (unlist (loghdr) < 0))
     stop ("error reading log: R does not support unsigned long integers.",
           "Please contact the maintainer of the package.")
   
	log <- list ()
	data <- list ()
	
	## read binary part of log
	if (log.bin) 
		log$.log.bin <- raw.data [loghdr$.last.read + seq_len (loghdr$logbins)]
	
	## read binary on-disk-only part of log
	if (log.disk) 
		log$.log.disk <- raw.data [loghdr$.last.read + loghdr$logbins + seq_len (loghdr$logdsks)]
	
	## read text part of log
	if (log.txt & loghdr$logsizd > loghdr$logtxto) {
		log.txt <- raw.data [pos + loghdr$logtxto + seq_len (loghdr$logsizd - loghdr$logtxto)]
      if (tail (log.txt, 1) ==  .nul)   # throw away nul at the end
        log.txt <- head (log.txt, -1)
      log.txt [log.txt == .nul] <- replace.nul
      log.txt <- readChar (log.txt, length (log.txt), useBytes=T)
      log.txt <- gsub (rawToChar (replace.nul), '\r\n', log.txt)
      log.txt <- iconv (log.txt, iconv.from, iconv.to)
#		log.txt <- paste (rawToChar (log.txt, multiple = TRUE), collapse = "")
		log.txt <- split.string (log.txt, "\r\n") ## spc file spec says \r\n regardless of OS
		log.txt <- split.line (log.txt, "=")
		## lapply (keys, function (x)
      ##         gsub (sprintf ("^.*\r\n%s\\s*=([^\r\n]+)\r\n.*$", x), "\\1", log.txt))
		data <- getbynames (log.txt, keys.log2data)
		log <- c (log, getbynames (log.txt, keys.log2log))
	}
	
	list (log.long = log, extra.data = data)
}


## read y data ......................................................................................
##

.spc.read.y <- function (raw.data, pos, npts, exponent, word) {
	if (exponent == -128) { # 4 byte float
		
		list (y = readBin (raw.data [pos + seq_len (npts * 4)], "numeric", npts, 4),
				.last.read = pos + npts * 4)
		
	} else if (word) { # 2 byte fixed point integer = word
		
		list (y = readBin (raw.data [pos + seq_len (npts * 2)], "integer", npts, 2, signed = TRUE) *
						2 ^ (exponent - 16),
				.last.read = pos + npts * 2)
		
	} else { # 4 byte fixed point integer = dword
		list (y = readBin (raw.data [pos + seq_len (npts * 4)], "integer", npts, 4) *
						2 ^ (exponent - 32),
				.last.read = pos + npts * 4)
	}
}

## read x data ......................................................................................
##

.spc.read.x <- function (raw.data, pos, npts) {
	list (x = readBin (raw.data [pos + seq_len (npts * 4)], "numeric", npts, 4),
			.last.read = pos + 4 * npts)
}

## error .............................................................................................
.spc.error <- function (fname, objects, ...) {
	cat ('ERROR in read.spc function ', fname, '\n\n')
	for (i in seq_along (objects)) {
		cat (names (objects) [i], ":\n")
		str (objects [[i]], vec.len = 20)
	}
	stop (...)
}

.spc.ftflags <- function (x) {
	ftflgs <- as.logical (x %/% 2^(0 : 7) %% 2)
	names (ftflgs) <- c ('TSPREC', 'TCGRAM', 'TMULTI', 'TRANDM',
			'TORDRD', 'TALABS', 'TXYXYS', 'TXVALS')
	ftflgs
}

#####################################################################################################


##' Import for Thermo Galactic's spc file format
##' These functions allow to import Thermo Galactic/Grams .spc files.
##' 
##' @param filename The complete file name of the .spc file.
##' @param keys.hdr2data,keys.log2data character vectors with the names of parameters in the .spc
##' file's log block (log2xxx) or header (hdr2xxx) that should go into the extra data (yyy2data) of
##' the returned hyperSpec object.
##'  
##' All header fields specified in the .spc file format specification (see
##'   below) are imported and can be referred to by their de-capitalized names.
##' @param keys.hdr2log,keys.log2log deprecated
##' @param log.txt Should the text part of the .spc file's log block be read?
##' @param log.bin,log.disk Should the normal and on-disk binary parts of the
##'   .spc file's log block be read?  If so, they will be put as raw vectors
##'   into the hyperSpec object's log.
##' @param hdr A list with fileheader fields that overwrite the settings of
##'   actual file's header.
##' 
##' Use with care, and look into the source code for detailed insight on the
##'   elements of this list.
##' @param no.object If \code{TRUE}, a list with wavelengths, spectra, labels,
##'   log and data are returned instead of a hyperSpec object.
##' 
##' This parameter will likely be subject to change in future - use with care.
##' @return If the file contains multiple spectra with individual wavelength
##'   axes, \code{read.spc} returns a list of hyperSpec objects.  Otherwise the
##'   result is a hyperSpec object.
##' 
##' \code{read.spc.KaiserMap} returns a hyperSpec object with data columns x,
##'   y, and z containing the stage position as recorded in the .spc files'
##'   log.
##' @note Only a restricted set of test files was available for development.
##'   Particularly, the w-planes feature could not be tested.
##' 
##' If you have .spc files that cannot be read with these function, don't
##'   hesitate to contact the package maintainer with your code patch or asking
##'   advice.
##' @author C. Beleites
##' @rdname read-spc
##' @seealso \code{\link[hyperSpec]{textio}}
##' @references Source development kit and file format specification of .spc
##'   files.
##' @export
##' @keywords IO file
##' @examples
##' 
##' ## get the sample .spc files from ftirsearch.com (see above)
##' \dontrun{
##' # single spectrum
##' spc <- read.spc ("BENZENE.SPC")
##' plot (spc)
##' 
##' # multi-spectra .spc file with common wavelength axis
##' spc <- read.spc ('IG_MULTI.SPC')
##' spc
##' 
##' # multi-spectra .spc file with individual wavelength axes
##' spc <- read.spc ("BARBITUATES.SPC")
##' plot (spc [[1]], lines.args = list (type = "h"))
##' }
##' 
read.spc <- function (filename,
		keys.hdr2data = FALSE, keys.hdr2log = FALSE,
		keys.log2data = FALSE, keys.log2log = FALSE,
		log.txt = TRUE, log.bin = FALSE, log.disk = FALSE,
		hdr = list (),
		no.object = FALSE){
	
	## f contains the raw bytes of the file
	
	## fpos marks the position of the last read byte
	## this is the same as the offset from beginning of the file (count 0) in the .spc definition  	
	
	f <- readBin (filename, "raw", file.info (filename)$size, 1)
	
	hdr <- modifyList (.spc.filehdr (f), hdr)
	fpos <- hdr$.last.read
	
	if (! hdr$ftflgs ['TXYXYS']) {
		if (! hdr$ftflgs ['TXVALS']) {
			## spectra with common evenly spaced wavelength axis
			wavelength <- seq (hdr$ffirst, hdr$flast, length.out = hdr$fnpts)
		} else {
			## spectra with common unevenly spaced wavelength axis
		#	if (! hdr$ftflgs ['TMULTI']) { # also for multifile with common wavelength axis
				tmp <- .spc.read.x (f, fpos, hdr$fnpts)
				wavelength <- tmp$x
				fpos <- tmp$.last.read
			}
		#}
	}
	
	## otherwise (TXYXYS set) hdr$fnpts gives offset to subfile directory if that exists
	
	## obtain labels from file hdr or from parameter	
	if (hdr$fwplanes > 0)
		label <- modifyList (list (w = hdr$fwtype), label)
	
	label <- list (.wavelength = hdr$fxtype, spc = hdr$fytype,
					z = hdr$fztype, z.end = hdr$fztype)
				
	## prepare list for hyperSpec log and data.frame for extra data
	
	data <- list (z = NA, z.end = NA)
	if (hdr$fwplanes > 0)
		data <- c (data, w = NA)
	
	## process the log block
	tmp <- .spc.log (f, hdr$flogoff,
			log.bin, log.disk, log.txt,
			keys.log2data,  keys.log2log)
	## TODO: remove keys.log2log data2log

	data <- c (data, tmp$extra.data, getbynames (hdr, keys.hdr2data))
	
	## try to preallocate spectra matrix and extra data data.frame
	## if multispectra file with separate wavelength axes, prepare a list
	if (hdr$ftflgs ['TXYXYS'] && hdr$ftflgs ['TMULTI']) {
		spc <- list ()		
		data <- as.data.frame (data, stringsAsFactors = FALSE)
	} else {
		spc <- matrix (NA, nrow = hdr$fnsub, ncol = hdr$fnpts)

      ## the *type header elements are expressions. They need to be converted to character.
      data <- sapply (data, function (x) {
        if (mode (x) == "expression")
            as.character (x)
        else
            x
      })
                       
      
		data <- as.data.frame (lapply (data, rep, hdr$fnsub), stringsAsFactors = FALSE)
	}
	
	## read subfiles
	if (hdr$subfiledir){ ## TXYXYS
		hdr$subfiledir <- .spc.subfiledir (f, hdr$subfiledir, hdr$fnsub)
		
		for (s in seq_len (hdr$fnsub)) {

			hdr <- .spc.subhdr (f, hdr$subfiledir$ssfposn [s], hdr)
			fpos <- hdr$.last.read
			wavelength <- .spc.read.x (f, fpos, hdr$subhdr$subnpts)
			fpos <- wavelength$.last.read
			
			y <- .spc.read.y (f, fpos, npts = hdr$subhdr$subnpts, exponent = hdr$subhdr$subexp,
					word = hdr$ftflgs ['TSPREC'])
			fpos <- y$.last.read
			
			data$z <- hdr$subhdr$subtime
			data$z.end <- hdr$subhdr$subnext
			
			if (hdr$fwplanes > 0)
				data$w <- hdr$subhdr$w
			
			if (! exists ('wavelength'))
				.spc.error ("read.spc", list (hdr = hdr),
						"wavelength not read. This may be caused by wrong header information.")
			
			spc [[s]] <- new ("hyperSpec",
					spc = y$y, 
					wavelength = wavelength$x,
					data = data,
					labels = label)
		}
		
	} else { ## multiple y data blocks behind each other
		for (s in seq_len (hdr$fnsub)) {
			hdr <- .spc.subhdr (f, fpos, hdr)
			fpos <- hdr$.last.read
			tmp <- .spc.read.y (f, fpos, npts = hdr$subhdr$subnpts, exponent = hdr$subhdr$subexp,
					word = hdr$ftflgs ['TSPREC'])
			fpos <- tmp$.last.read
			
			spc [s, ] <- tmp$y
			data [s, c('z', 'z.end')] <- unlist (hdr$subhdr [c('subtime', 'subnext')])
			
			if (hdr$fwplanes > 0)
				data [s, "w"] <- hdr$subhdr$w
		}
	}
	
	if (hdr$ftflgs ['TXYXYS'] && hdr$ftflgs ['TMULTI']) 
		spc
	else if (no.object)
		list (spc = spc, wavelength = wavelength, data = data, labels = label)
	else 
		new ("hyperSpec",  spc = spc, wavelength = wavelength,
				data = data [rep (1, hdr$fnsub), ], labels = label)
}


