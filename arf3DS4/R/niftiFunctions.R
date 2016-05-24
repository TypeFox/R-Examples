#############################################
# arf3DS4 S4 NIFTI FUNCTIONS				#
# Wouter D. Weeda							#
# University of Amsterdam					#
#############################################

#[CONTAINS]
#getFileInfo
#newFile
#readHeader
#readData				[user]
#writeHeader
#writeData				[user]
#writeHeaderPart
#writeDataPart
#headToName
#getIntent
#fmri2array				[user]
#bin2dec
#xyzt2char


getFileInfo <- 
function(filename)
# getFileInfo sets file information (file path, location, extension, gzippedness and endianness)
{
	#set separator
	sp <- .Platform$file.sep
	
	#check if only file else pre-append working directory
	if(length(grep(sp,filename))==0) filename <- paste(getwd(),sp,filename,sep='')
	
	#check if file exists
	if(!file.exists(filename)) stop(paste('[getFileInfo] File',filename,'does noet exist!'))
	
	#create new instance of class fileinfo
	fileinf=new('nifti.fileinfo')
	
	#split filename into path, name and extensions
	fn <- strsplit(filename,sp)[[1]][length(strsplit(filename,sp)[[1]])]
	#fsplit_name <- strsplit(tolower(fn),'\\.')[[1]]
	fsplit_name <- strsplit(fn,'\\.')[[1]]
	fsplit_path <- gsub(paste(sp,fn,sep=''),'',filename)
	fname <- strsplit(fsplit_name[1],.Platform$file.sep)[[1]]	
	
	#set filename and path
	.nifti.fileinfo.filename(fileinf) <- fname[length(fname)]
	.nifti.fileinfo.fullpath(fileinf) <- fsplit_path
	
	#check gzippedness and set extension
	if(fsplit_name[length(fsplit_name)]=='gz') {
		.nifti.fileinfo.gzipped(fileinf) <- TRUE
		.nifti.fileinfo.extension(fileinf) <- fsplit_name[length(fsplit_name)-1]	
	} else {
		.nifti.fileinfo.gzipped(fileinf) <- FALSE
		.nifti.fileinfo.extension(fileinf) <- fsplit_name[length(fsplit_name)]
	}
	
	#check valid files if img/hdr pair and if IMG try and open header	
	if(.nifti.fileinfo.extension(fileinf)=='img') {
		if(.nifti.fileinfo.gzipped(fileinf)==TRUE) {
			if(file.exists(paste(.nifti.fileinfo.fullpath(fileinf),sp,.nifti.fileinfo.filename(fileinf),'.hdr.gz',sep=''))) {
				.nifti.fileinfo.extension(fileinf) <- 'hdr'	
			} else {
				stop('No valid img/hdr pair found. HDR does not exist.\n')
			}
		} else {
			if(file.exists(paste(.nifti.fileinfo.fullpath(fileinf),sp,.nifti.fileinfo.filename(fileinf),'.hdr',sep=''))) {
				.nifti.fileinfo.extension(fileinf) <- 'hdr'	
			} else {
				stop('No valid img/hdr pair found. HDR does not exist.\n')
			}
		} 
	}
	
	if(.nifti.fileinfo.extension(fileinf)=='hdr') {
		if(.nifti.fileinfo.gzipped(fileinf)==TRUE) {
			if(!file.exists(paste(.nifti.fileinfo.fullpath(fileinf),sp,.nifti.fileinfo.filename(fileinf),'.img.gz',sep=''))) stop('No valid img/hdr pair found. IMG does not exist.\n')	
		} else {
			if(!file.exists(paste(.nifti.fileinfo.fullpath(fileinf),sp,.nifti.fileinfo.filename(fileinf),'.img',sep=''))) stop('No valid img/hdr pair found. IMG does noet exist.\n')	
		}
	}
	
	#check nifti header. Returns fileinf of class nifti.header
	headinfo=readHeader(fileinf)
	
	#check if header size is 348 (if not try other endian)
	if(.nifti.header.sizeof_hdr(headinfo)!=348) {
		
		#change endian
		if(.nifti.fileinfo.endian(fileinf)=='big') .nifti.fileinfo.endian(fileinf) <- 'little' else .nifti.fileinfo.endian(fileinf) <- 'big'
		
		#check header again
		headinfo=readHeader(fileinf)
		
		if(.nifti.header.sizeof_hdr(headinfo)!=348) {
			#header is of incorrect size
			stop('Header is not of size 348 (even after swapping endian). File might not be of the correct format.\n')
		}
	}
	
	#set filetype according to header info (not based on extension)
	if(.nifti.header.magic(headinfo)=='n+1') {
		.nifti.header.filetype(headinfo) <- 'nifti+1'
	} else	if(.nifti.header.magic(headinfo)=='ni1') {
		.nifti.header.filetype(headinfo) <- 'nifti1'
		
	} else if(.nifti.header.magic(headinfo)=='') {
		.nifti.header.filetype(headinfo) <- 'analyze'
		
	} else {
		stop('Magicstring contains unknown characters. File might not be the correct format.\n')
	}
	
	#Returns object of class nifti.header		
	return(invisible(headinfo))
	
}

newFile <- 
function(filename,templateHDR) 
#set header info accordomg to a templateHeader
{
	##newFile creates a new file of the given filename (and directories it should be in if they not exist)
	## input is a filename (full) and template headerinfo
	## output is headerInfo of the newly created file
	
	sp <- .Platform$file.sep
	
	#split filename into path, name and extensions
	fn <- strsplit(filename,sp)[[1]][length(strsplit(filename,sp)[[1]])]
	#fsplit_name <- strsplit(tolower(fn),'\\.')[[1]]
	fsplit_name <- strsplit(fn,'\\.')[[1]]
	fsplit_path <- gsub(paste(sp,fn,sep=''),'',filename)
	fname <- strsplit(fsplit_name[1],.Platform$file.sep)[[1]]	
	
	#set filename and path
	.nifti.header.filename(templateHDR) <- fname[length(fname)]
	.nifti.header.fullpath(templateHDR) <- fsplit_path
	
	#check gzippedness and set extension
	if(fsplit_name[length(fsplit_name)]=='gz') {
		.nifti.header.gzipped(templateHDR) <- TRUE
		.nifti.header.extension(templateHDR) <- fsplit_name[length(fsplit_name)-1]	
	} else {
		.nifti.header.gzipped(templateHDR) <- FALSE
		.nifti.header.extension(templateHDR) <- fsplit_name[length(fsplit_name)]
	}
	
	# return header info object for created file
	return(templateHDR)
}

readHeader <- 
function(fileinf) 
# readHeader reads nifti header info from nii,img,hdr and gzipped exts.
{
	sp <- .Platform$file.sep
	
	#open new header headinf
	headinf=new('nifti.header',fileinf)
	
	#open based on gzippedness
	if(.nifti.header.gzipped(headinf)==TRUE) {
		con <- gzfile(paste(.nifti.header.fullpath(headinf),sp,.nifti.header.filename(headinf),'.',.nifti.header.extension(headinf),'.gz',sep=''),open='rb')	
	} else {
		con <- file(paste(.nifti.header.fullpath(headinf),sp,.nifti.header.filename(headinf),'.',.nifti.header.extension(headinf),sep=''),open='rb')		
	}
	
	#if connection is open try and read header info (magicstring and headersize)
	if(con) {
		
		readin <- function(con,headinf) {
			#read in all elements (independent of NIFTI or ANALYZE), strings with nulls are truncated	  	  							 #byte offset, type, descr.
			.nifti.header.sizeof_hdr(headinf) <- readBin(con,integer(),n=1,size=4,signed=T,endian=.nifti.header.endian(headinf))         #  0, int32,  MUST BE 348
			.nifti.header.data_type(headinf) <- rawToChar(removeNullString(readBin(con,raw(),n=10)))                                	 #  4, char[10]
			.nifti.header.db_name(headinf) <- rawToChar(removeNullString(readBin(con,raw(),n=18)))                                  	 # 14, char[18]
			.nifti.header.extents(headinf) <- readBin(con,integer(),n=1,size=4,signed=T,endian=.nifti.header.endian(headinf))            # 32, int32
			.nifti.header.session_error(headinf) <- readBin(con,integer(),n=1,size=2,signed=T,endian =.nifti.header.endian(headinf))     # 36, short
			.nifti.header.regular(headinf) <- rawToChar(removeNullString(readBin(con,raw(),n=1)))                                   	 # 38, char[1]
			.nifti.header.dim_info(headinf) <- rawToChar(removeNullString(readBin(con,raw(),n=1)))                                  	 # 39, char[1], MRI SLICE ORDERING
			.nifti.header.dims(headinf) <- readBin(con,integer(),n=8,size=2,signed=T,endian=.nifti.header.endian(headinf))               # 40, short[8], DATA ARRAY DIMENSIONS
			.nifti.header.intent_p1(headinf) <- readBin(con,double(),n=1,size=4,endian=.nifti.header.endian(headinf))                    # 56, float, 1ST INTENT PARAMETER
			.nifti.header.intent_p2(headinf) <- readBin(con,double(),n=1,size=4,endian=.nifti.header.endian(headinf))                    # 60, float, 2ND INTENT PARAMETER
			.nifti.header.intent_p3(headinf) <- readBin(con,double(),n=1,size=4,endian=.nifti.header.endian(headinf))                    # 64, float, 3RD INTENT PARAMETER
			.nifti.header.intent_code(headinf) <- readBin(con,integer(),n=1,size=2,signed=T,endian=.nifti.header.endian(headinf))        # 68, short, NIFITINTENT CODE
			.nifti.header.datatype(headinf) <- readBin(con,integer(),n=1,size=2,signed=T,endian=.nifti.header.endian(headinf))           # 70, short, DEFINES DATA TYPE
			.nifti.header.bitpix(headinf) <- readBin(con,integer(),n=1,size=2,signed=T,endian=.nifti.header.endian(headinf))             # 72, short, NUMBER BITS PER VOXEL
			.nifti.header.slice_start(headinf) <- readBin(con,integer(),n=1,size=2,signed=T,endian=.nifti.header.endian(headinf))        # 74, short, FIRST SLICE INDEX
			.nifti.header.pixdim(headinf) <- readBin(con,double(),n=8,size=4,endian=.nifti.header.endian(headinf))                       # 76, float[8], GRID SPACINGS
			.nifti.header.vox_offset(headinf) <- readBin(con,double(),n=1,size=4,endian=.nifti.header.endian(headinf))                   #108, float, OFFSET INTO .NII FILE
			.nifti.header.scl_slope(headinf) <- readBin(con,double(),n=1,size=4,endian=.nifti.header.endian(headinf))                    #112, float, DATA SCALING: SLOPE
			.nifti.header.scl_inter(headinf) <- readBin(con,double(),n=1,size=4,endian=.nifti.header.endian(headinf))                    #116, float, DATA SCALING: INTERCEPT
			.nifti.header.slice_end (headinf) <- readBin(con,integer(),n=1,size=2,signed=T,endian=.nifti.header.endian(headinf))         #120, short, LAST SLICE INDEX
			.nifti.header.slice_code(headinf) <- readBin(con,integer(),size=1,n=1,signed=T,endian=.nifti.header.endian(headinf))         #122, integer[1], SLICE TIMING ORDER
			#.nifti.header.xyzt_units(headinf) <- readBin(con,integer(),size=1,n=1,signed=T,endian=.nifti.header.endian(headinf))        #123, raw[1], UNITS OF PIXDIM
			.nifti.header.xyzt_units(headinf) <- rawToChar(removeNullString(readBin(con,raw(),n=1)))									 #123, raw[1], UNITS OF PIXDIM	
			.nifti.header.cal_max(headinf) <- readBin(con,double(),n=1,size=4,endian=.nifti.header.endian(headinf))                      #124, float, MAX DISPLAY INTENSITY
			.nifti.header.cal_min(headinf) <- readBin(con,double(),n=1,size=4,endian=.nifti.header.endian(headinf))                      #128, float, MIN DISPLAY INTENSITY
			.nifti.header.slice_duration(headinf) <- readBin(con,double(),n=1,size=4,endian=.nifti.header.endian(headinf))               #132, float, TIME FOR 1 SLICE
			.nifti.header.toffset(headinf) <- readBin(con,double(),n=1,size=4,endian=.nifti.header.endian(headinf))                      #136, float, TIME AXIS SHIFT
			.nifti.header.glmax(headinf) <- readBin(con,integer(),n=1,size=4,signed=T,endian=.nifti.header.endian(headinf))              #140, int32
			.nifti.header.glmin(headinf) <- readBin(con,integer(),n=1,size=4,signed=T,endian=.nifti.header.endian(headinf))              #144, int32
			.nifti.header.descrip(headinf) <- rawToChar(removeNullString(readBin(con,raw(),n=80))) 		  								 #148, char[80], TEXT
			.nifti.header.aux_file(headinf) <- rawToChar(removeNullString(readBin(con,raw(),n=24)))		  								 #228, char[24], AUXILIARY FILENAME
			.nifti.header.qform_code(headinf) <- readBin(con,integer(),n=1,size=2,signed=T,endian=.nifti.header.endian(headinf))         #252, short, NIFITXFORM CODE
			.nifti.header.sform_code(headinf) <- readBin(con,integer(),n=1,size=2,signed=T,endian=.nifti.header.endian(headinf))         #254, short, NIFITXFORM CODE
			.nifti.header.quatern_b(headinf) <- readBin(con,double(),n=1,size=4,endian=.nifti.header.endian(headinf))                    #256, float, QUATERNION B PARAM
			.nifti.header.quatern_c(headinf) <- readBin(con,double(),n=1,size=4,endian=.nifti.header.endian(headinf))                    #260, float, QUATERNION C PARAM
			.nifti.header.quatern_d(headinf) <- readBin(con,double(),n=1,size=4,endian=.nifti.header.endian(headinf))                    #264, float, QUATERNION D PARAM
			.nifti.header.qoffset_x(headinf) <- readBin(con,double(),n=1,size=4,endian=.nifti.header.endian(headinf))                    #268, float, QUATERNION X SHIFT
			.nifti.header.qoffset_y(headinf) <- readBin(con,double(),n=1,size=4,endian=.nifti.header.endian(headinf))                    #272, float, QUATERNION Y SHIFT
			.nifti.header.qoffset_z(headinf) <- readBin(con,double(),n=1,size=4,endian=.nifti.header.endian(headinf))                    #276, float, QUATERNION Z SHIFT
			.nifti.header.srow_x(headinf) <- readBin(con,double(),n=4,size=4,endian=.nifti.header.endian(headinf))                       #280, float[4], 1ST ROW AFFINE TRANSFORM
			.nifti.header.srow_y(headinf) <- readBin(con,double(),n=4,size=4,endian=.nifti.header.endian(headinf))                       #296, float[4], 2ND ROW AFFINE TRANSFORM
			.nifti.header.srow_z(headinf) <- readBin(con,double(),n=4,size=4,endian=.nifti.header.endian(headinf))                       #312, float[4], 3RD ROW AFFINE TRANSFORM
			.nifti.header.intent_name(headinf) <- rawToChar(removeNullString(readBin(con,raw(),n=16)))                              		  								#328, char[16], NAME OR MEANING OF DATA
			.nifti.header.magic(headinf) <- rawToChar(removeNullString(readBin(con,raw(),n=4)))    																		#344, char[4], MAGICSTRING!
		
			#close connection
			close(con)

			#set fmri data attributes for read in
			if(.nifti.header.datatype(headinf)==0) { .nifti.header.data.type(headinf)<-'none';.nifti.header.data.signed(headinf)<-TRUE } else
			if(.nifti.header.datatype(headinf)==1) { .nifti.header.data.type(headinf)<-'raw';.nifti.header.data.signed(headinf)<-TRUE } else
			if(.nifti.header.datatype(headinf)==2) { .nifti.header.data.type(headinf)<-'integer';.nifti.header.data.signed(headinf)<-FALSE } else
			if(.nifti.header.datatype(headinf)==4) { .nifti.header.data.type(headinf)<-'integer';.nifti.header.data.signed(headinf)<-TRUE } else
			if(.nifti.header.datatype(headinf)==8) { .nifti.header.data.type(headinf)<-'integer';.nifti.header.data.signed(headinf)<-TRUE } else
			if(.nifti.header.datatype(headinf)==16) { .nifti.header.data.type(headinf)<-'double';.nifti.header.data.signed(headinf)<-TRUE } else
			if(.nifti.header.datatype(headinf)==32) { .nifti.header.data.type(headinf)<-'complex';.nifti.header.data.signed(headinf)<-TRUE } else
			if(.nifti.header.datatype(headinf)==64) { .nifti.header.data.type(headinf)<-'double';.nifti.header.data.signed(headinf)<-TRUE } else
			if(.nifti.header.datatype(headinf)==256) { .nifti.header.data.type(headinf)<-'integer';.nifti.header.data.signed(headinf)<-TRUE } else
			if(.nifti.header.datatype(headinf)==512) { .nifti.header.data.type(headinf)<-'integer';.nifti.header.data.signed(headinf)<-FALSE } else
			if(.nifti.header.datatype(headinf)==768) { .nifti.header.data.type(headinf)<-'integer';.nifti.header.data.signed(headinf)<-FALSE } else
			if(.nifti.header.datatype(headinf)==1024) { .nifti.header.data.type(headinf)<-'integer';.nifti.header.data.signed(headinf)<-TRUE } else
			if(.nifti.header.datatype(headinf)==1280) { .nifti.header.data.type(headinf)<-'integer';.nifti.header.data.signed(headinf)<-FALSE } else
			if(.nifti.header.datatype(headinf)==1536) { .nifti.header.data.type(headinf)<-'double';.nifti.header.data.signed(headinf)<-TRUE } else
			if(.nifti.header.datatype(headinf)==1792) { .nifti.header.data.type(headinf)<-'complex';.nifti.header.data.signed(headinf)<-TRUE } else
			if(.nifti.header.datatype(headinf)==2048) { .nifti.header.data.type(headinf)<-'complex';.nifti.header.data.signed(headinf)<-TRUE } else { stop('Datatype is unknown!') }
			
			return(headinf)
		}
		
		#Read in Data with warning suppression	
		headinf <- suppressWarnings(readin(con,headinf))
		
	} else {
		stop('Unable to open connection.\n')
	}

	
	#return nifti.header object
	return(invisible(headinf))

	
}

readData <-
function(filename) 
# readData reads in nifti/analyze datafiles
{
	sp <- .Platform$file.sep
	
	#obtain header info
	headinf <- readHeader(getFileInfo(filename))
	
	## set correct extension
	if(.nifti.header.extension(headinf)=='hdr') extension='img' else extension <- .nifti.header.extension(headinf)
	
	#open based on gzippedness
	if(.nifti.header.gzipped(headinf)==TRUE) {
		con <- gzfile(paste(.nifti.header.fullpath(headinf),sp,.nifti.header.filename(headinf),'.',extension,'.gz',sep=''),open='rb')	
	} else {
		con <- file(paste(.nifti.header.fullpath(headinf),sp,.nifti.header.filename(headinf),'.',extension,sep=''),open='rb')		
	}
	
	if(con) {
		#create new vecdat object
		data <- new('fmri.data',headinf)
		
		#set length of data to read in
		n=(.nifti.header.dims(headinf)[2]*.nifti.header.dims(headinf)[3]*.nifti.header.dims(headinf)[4]*.nifti.header.dims(headinf)[5]*(.nifti.header.bitpix(headinf)/8))
		
		#read everthing before vox_offset
		readBin(con,raw(),.nifti.header.vox_offset(headinf))
		
		#read in data
		.fmri.data.datavec(data) <- readBin(con, what=.nifti.header.data.type(headinf), n=n, size=(.nifti.header.bitpix(headinf)/8), signed=.nifti.header.data.signed(headinf), endian=.nifti.header.endian(headinf))
		
	} else stop('Unable to open connection') 	
	
	#close connections
	close(con)
	
	#return datavector
	return(invisible(data))

}

writeHeader <- 
function(headinf) 
#writeHeader writes header info to a .hdr file
{
	sp <- .Platform$file.sep
	
	.nifti.header.vox_offset(headinf) <- 0
	
	## set correct extension
	if(.nifti.header.extension(headinf)=='img') extension='hdr' else extension <- .nifti.header.extension(headinf)
	
	## open based on gzippedness
	if(.nifti.header.gzipped(headinf)==TRUE) {
		fn <- paste(.nifti.header.fullpath(headinf),sp,.nifti.header.filename(headinf),'.',extension,'.gz',sep='')
		con <- gzfile(fn,open='wb')	
	} else {
		fn <-paste(.nifti.header.fullpath(headinf),sp,.nifti.header.filename(headinf),'.',extension,sep='')
		con <- file(fn,open='wb')	
	}
	
	## write to .hdr file
	if(extension=='hdr') con <- writeHeaderPart(con,headinf) else warning('Header info can only be written to .hdr files! No header info is written!')
	
	close(con)
	
	return(invisible(TRUE))
	
}

writeData <- 
function(headinf,datavec=NULL) 
#write data to a nifti file
{
	## writeData writes both Header and Data info to a file
	## input is headerinfo and datavector
	## output is logical
	sp <- .Platform$file.sep
	
	if(is.null(datavec)) datavec = .fmri.data.datavec(headinf)
	
	# set correct extension
	if(.nifti.header.extension(headinf)=='hdr') extension='img' else extension <- .nifti.header.extension(headinf)
	
	#open files based on gzippedness 
	if(.nifti.header.gzipped(headinf)==TRUE) {
		fn <- paste(.nifti.header.fullpath(headinf),sp,.nifti.header.filename(headinf),'.',extension,'.gz',sep='')
		con <- gzfile(fn,open='wb')	
	} else {
		fn <-paste(.nifti.header.fullpath(headinf),sp,.nifti.header.filename(headinf),'.',extension,sep='')
		con <- file(fn,open='wb')	
	}
	
	#set correct cal_min and cal_max
	.nifti.header.cal_max(headinf) = max(datavec)
	.nifti.header.cal_min(headinf) = min(datavec)
		
	# for nii files write header and data in one connection (set offset of data as 352)
	if(extension=='nii') {
		.nifti.header.vox_offset(headinf) <- 352
		con <- writeHeaderPart(con,headinf)
		con <- writeDataPart(con,headinf,datavec)
		
	}
	
	# for img/hdr pairs write header to .hdr and data to .img
	if(extension=='img') {
		.nifti.header.vox_offset(headinf) <- 0
		writeHeader(headinf)
		con <- writeDataPart(con,headinf,datavec)
	}
	
	#close connection
	close(con)
	
	return(TRUE)
	
}

writeHeaderPart <- 
function(con,headinf) 
#writeHeaderPart of a nifti file
{
	## writeHeaderPart writes nifti/analyze headerfiles
	## input is connection,headerinfo,and datavec
	## output is connection
	
	#if connection is open try and read header info (magicstring and headersize)
	if(con) {
		#write all elements (independent of NIFTI or ANALYZE), watch String Truncation		  		  	  			   #byte offset, type, descr.
		writeBin(as.integer(.nifti.header.sizeof_hdr(headinf)),con,size=4,endian=.nifti.header.endian(headinf))		   #  0, int32,  MUST BE 348
		writeBin(charToRaw(.nifti.header.data_type(headinf)),con,endian=.nifti.header.endian(headinf))                 #  4, char[10]
		writeBin(raw(10-nchar(.nifti.header.data_type(headinf))),con,endian=.nifti.header.endian(headinf)) 			
		writeBin(charToRaw(.nifti.header.db_name(headinf)),con,endian=.nifti.header.endian(headinf))                   # 14, char[18]
		writeBin(raw(18-nchar(.nifti.header.db_name(headinf))),con,endian=.nifti.header.endian(headinf)) 
		writeBin(as.integer(.nifti.header.extents(headinf)),con,size=4,endian=.nifti.header.endian(headinf))           # 32, int32
		writeBin(as.integer(.nifti.header.session_error(headinf)),con,size=2,endian=.nifti.header.endian(headinf))     # 36, short
		writeBin(charToRaw(.nifti.header.regular(headinf)),con,endian=.nifti.header.endian(headinf))                   # 38, char[1]
		writeBin(raw(1-nchar(.nifti.header.regular(headinf))),con,endian=.nifti.header.endian(headinf)) 
		writeBin(charToRaw(.nifti.header.dim_info(headinf)),con,endian=.nifti.header.endian(headinf))                  # 39, char[1], MRI SLICE ORDERING
		writeBin(raw(1-nchar(.nifti.header.dim_info(headinf))),con,endian=.nifti.header.endian(headinf)) 
		writeBin(as.integer(.nifti.header.dims(headinf)),con,size=2,endian=.nifti.header.endian(headinf))              # 40, short[8], DATA ARRAY DIMENSIONS
		writeBin(as.double(.nifti.header.intent_p1(headinf)),con,size=4,endian=.nifti.header.endian(headinf))          # 56, float, 1ST INTENT PARAMETER
		writeBin(as.double(.nifti.header.intent_p2(headinf)),con,size=4,endian=.nifti.header.endian(headinf))          # 60, float, 2ND INTENT PARAMETER
		writeBin(as.double(.nifti.header.intent_p3(headinf)),con,size=4,endian=.nifti.header.endian(headinf))          # 64, float, 3RD INTENT PARAMETER
		writeBin(as.integer(.nifti.header.intent_code(headinf)),con,size=2,endian=.nifti.header.endian(headinf))       # 68, short, NIFITINTENT CODE
		writeBin(as.integer(.nifti.header.datatype(headinf)),con,size=2,endian=.nifti.header.endian(headinf))          # 70, short, DEFINES DATA TYPE
		writeBin(as.integer(.nifti.header.bitpix(headinf)),con,size=2,endian=.nifti.header.endian(headinf))            # 72, short, NUMBER BITS PER VOXEL
		writeBin(as.integer(.nifti.header.slice_start(headinf)),con,size=2,endian=.nifti.header.endian(headinf))       # 74, short, FIRST SLICE INDEX
		writeBin(as.double(.nifti.header.pixdim(headinf)),con,size=4,endian=.nifti.header.endian(headinf))             # 76, float[8], GRID SPACINGS
		writeBin(as.double(.nifti.header.vox_offset(headinf)),con,size=4,endian=.nifti.header.endian(headinf))         #108, float, OFFSET INTO .NII FILE
		writeBin(as.double(.nifti.header.scl_slope(headinf)),con,size=4,endian=.nifti.header.endian(headinf))          #112, float, DATA SCALING: SLOPE
		writeBin(as.double(.nifti.header.scl_inter(headinf)),con,size=4,endian=.nifti.header.endian(headinf))          #116, float, DATA SCALING: INTERCEPT
		writeBin(as.integer(.nifti.header.slice_end (headinf)),con,size=2,endian=.nifti.header.endian(headinf))        #120, short, LAST SLICE INDEX
		writeBin(as.integer(.nifti.header.slice_code(headinf)),con,size=1,endian=.nifti.header.endian(headinf))        #122, integer[1], SLICE TIMING ORDER
		#writeBin(as.integer(.nifti.header.xyzt_units(headinf)),con,size=1,endian=.nifti.header.endian(headinf))        #123, integer[1], UNITS OF PIXDIM
		writeBin(charToRaw(.nifti.header.xyzt_units(headinf)),con,endian=.nifti.header.endian(headinf))       			#123, integer[1], UNITS OF PIXDIM	
		writeBin(as.double(.nifti.header.cal_max(headinf)),con,size=4,endian=.nifti.header.endian(headinf))            #124, float, MAX DISPLAY INTENSITY
		writeBin(as.double(.nifti.header.cal_min(headinf)),con,size=4,endian=.nifti.header.endian(headinf))            #128, float, MIN DISPLAY INTENSITY
		writeBin(as.double(.nifti.header.slice_duration(headinf)),con,size=4,endian=.nifti.header.endian(headinf))     #132, float, TIME FOR 1 SLICE
		writeBin(as.double(.nifti.header.toffset(headinf)),con,size=4,endian=.nifti.header.endian(headinf))            #136, float, TIME AXIS SHIFT
		writeBin(as.integer(.nifti.header.glmax(headinf)),con,size=4,endian=.nifti.header.endian(headinf))             #140, int32
		writeBin(as.integer(.nifti.header.glmin(headinf)),con,size=4,endian=.nifti.header.endian(headinf))             #144, int32
		writeBin(charToRaw(.nifti.header.descrip(headinf)),con,endian=.nifti.header.endian(headinf))                   #148, char[80], TEXT
		writeBin(raw(80-nchar(.nifti.header.descrip(headinf))),con,endian=.nifti.header.endian(headinf)) 
		writeBin(charToRaw(.nifti.header.aux_file(headinf)),con,endian=.nifti.header.endian(headinf))                  #228, char[24], AUXILIARY FILENAME
		writeBin(raw(24-nchar(.nifti.header.aux_file(headinf))),con,endian=.nifti.header.endian(headinf)) 
		writeBin(as.integer(.nifti.header.qform_code(headinf)),con,size=2,endian=.nifti.header.endian(headinf))        #252, short, NIFITXFORM CODE
		writeBin(as.integer(.nifti.header.sform_code(headinf)),con,size=2,endian=.nifti.header.endian(headinf))        #254, short, NIFITXFORM CODE
		writeBin(as.double(.nifti.header.quatern_b(headinf)),con,size=4,endian=.nifti.header.endian(headinf))          #256, float, QUATERNION B PARAM
		writeBin(as.double(.nifti.header.quatern_c(headinf)),con,size=4,endian=.nifti.header.endian(headinf))          #260, float, QUATERNION C PARAM
		writeBin(as.double(.nifti.header.quatern_d(headinf)),con,size=4,endian=.nifti.header.endian(headinf))          #264, float, QUATERNION D PARAM
		writeBin(as.double(.nifti.header.qoffset_x(headinf)),con,size=4,endian=.nifti.header.endian(headinf))          #268, float, QUATERNION X SHIFT
		writeBin(as.double(.nifti.header.qoffset_y(headinf)),con,size=4,endian=.nifti.header.endian(headinf))          #272, float, QUATERNION Y SHIFT
		writeBin(as.double(.nifti.header.qoffset_z(headinf)),con,size=4,endian=.nifti.header.endian(headinf))          #276, float, QUATERNION Z SHIFT
		writeBin(as.double(.nifti.header.srow_x(headinf)),con,size=4,endian=.nifti.header.endian(headinf))             #280, float[4], 1ST ROW AFFINE TRANSFORM
		writeBin(as.double(.nifti.header.srow_y(headinf)),con,size=4,endian=.nifti.header.endian(headinf))             #296, float[4], 2ND ROW AFFINE TRANSFORM
		writeBin(as.double(.nifti.header.srow_z(headinf)),con,size=4,endian=.nifti.header.endian(headinf))             #312, float[4], 3RD ROW AFFINE TRANSFORM
		writeBin(charToRaw(.nifti.header.intent_name(headinf)),con,endian=.nifti.header.endian(headinf))               #328, char[16], NAME OR MEANING OF DATA
		writeBin(raw(16-nchar(.nifti.header.intent_name(headinf))),con,endian=.nifti.header.endian(headinf)) 			
		writeBin(charToRaw(.nifti.header.magic(headinf)),con,endian=.nifti.header.endian(headinf))      			   #344, char[4], MAGICSTRING!
		writeBin(raw(4-nchar(.nifti.header.magic(headinf))),con,endian=.nifti.header.endian(headinf)) 
		
		if(.nifti.header.vox_offset(headinf)>348) for(i in 1:(.nifti.header.vox_offset(headinf)-348))  writeBin(as.integer(0),size=1,con,endian=.nifti.header.endian(headinf))
	
		
	} else {
		stop('Unable to open connection.\n')
	}
	
	#return nifti.header object
	return(con)
	
	
}

writeDataPart <- 
function(con,headinf,datavec) 
#write a datavector to nifti file
{
	## writeDataPart writes nifti/analyze datafiles
	## input is connection,headerinfo,and datavec
	## output is connection
			
	#write data to connection
	if(con) {
		mode(datavec) <- .nifti.header.data.type(headinf)	
		writeBin(datavec,con,size=(.nifti.header.bitpix(headinf)/8))
	} else stop('Unable to open connection') 

	
	#return TRUE
	return(con)

}

headToName <- 
function(headinf) 
#turn header information to a filename
{
	
	#if extension is nii
	if(.nifti.header.extension(headinf)=='nii') {
		if(.nifti.header.gzipped(headinf)) filename <- paste(.nifti.header.filename(headinf),'.nii.gz',sep='') else filename <- paste(.nifti.header.filename(headinf),'.nii',sep='')
	}
	
	#if extension is hdr
	if(.nifti.header.extension(headinf)=='hdr') {
		if(.nifti.header.gzipped(headinf)) filename <- paste(.nifti.header.filename(headinf),'.hdr.gz',sep='') else filename <- paste(.nifti.header.filename(headinf),'.hdr',sep='')
	}
	
	#if extension is img
	if(.nifti.header.extension(headinf)=='img') {
		if(.nifti.header.gzipped(headinf)) filename <- paste(.nifti.header.filename(headinf),'.img.gz',sep='') else filename <- paste(.nifti.header.filename(headinf),'.img',sep='')
	}
	
	#return filename
	return(filename)
		
}


getIntent <- 
function(code) 
#get Nifti Intent code
{
	intent = 'Unkown code'
	
	#statistics
	if(as.numeric(code)==2) intent = 'correlation'
	if(as.numeric(code)==3) intent = 't-test'
	if(as.numeric(code)==4) intent = 'F-test'
	if(as.numeric(code)==5) intent = 'z-score'
	if(as.numeric(code)==6) intent = 'chisq'
	if(as.numeric(code)==7) intent = 'beta'
	if(as.numeric(code)==8) intent = 'binom'
	if(as.numeric(code)==9) intent = 'gamma'
	if(as.numeric(code)==10) intent = 'poisson'
	if(as.numeric(code)==11) intent = 'normal'
	if(as.numeric(code)==12) intent = 'F-test nonc'
	if(as.numeric(code)==13) intent = 'chisq nonc'
	if(as.numeric(code)==14) intent = 'logistic'
	if(as.numeric(code)==15) intent = 'laplace'
	if(as.numeric(code)==16) intent = 'uniform'
	if(as.numeric(code)==17) intent = 't-test nonc'
	if(as.numeric(code)==18) intent = 'weibull'
	if(as.numeric(code)==19) intent = 'Chi'
	if(as.numeric(code)==20) intent = 'invgauss'
	if(as.numeric(code)==21) intent = 'extval (p-val)'
	if(as.numeric(code)==22) intent = 'p-value'
	if(as.numeric(code)==23) intent = 'logpval'
	if(as.numeric(code)==24) intent = 'log10pval'
	
	#Non-stat
	if(as.numeric(code)==1001) intent = 'estimate'
	if(as.numeric(code)==1002) intent = 'label'
	if(as.numeric(code)==1003) intent = 'neuroname'
	if(as.numeric(code)==1004) intent = 'genmatrix'
	if(as.numeric(code)==1005) intent = 'symmatrix'
	if(as.numeric(code)==1006) intent = 'displacement'
	if(as.numeric(code)==1007) intent = 'vector'
	if(as.numeric(code)==1008) intent = 'pointset'
	if(as.numeric(code)==1009) intent = 'triangle'
	if(as.numeric(code)==1010) intent = 'quaternion'
	if(as.numeric(code)==1011) intent = 'dimless'
	
	#gifti
	if(as.numeric(code)==2001) intent = 'timeseries'
	if(as.numeric(code)==2002) intent = 'node'
	if(as.numeric(code)==2003) intent = 'rgb-triplet'
	if(as.numeric(code)==2004) intent = 'rgba-vector'
	if(as.numeric(code)==2005) intent = 'shape'
	

	return(intent)
}


fmri2array <- function(fmridat) 
#converts a fmri.data object to an array
{
	dx = .fmri.data.dims(fmridat)[2]
	dy = .fmri.data.dims(fmridat)[3]
	dz = .fmri.data.dims(fmridat)[4]
	dt = .fmri.data.dims(fmridat)[5]
	
	out = .fmri.data.datavec(fmridat)
	
	if(.fmri.data.dims(fmridat)[1]==3) dim(out) = c(dx,dy,dz)
	if(.fmri.data.dims(fmridat)[1]==4) dim(out) = c(dx,dy,dz,dt)
	
	return(out)
}

bin2dec <- 
function(binvec) 
#calculate decimal value for a binary vector
{
	value = as.integer(binvec[1])
	if(length(binvec)>1) {
		
		for(i in 2:(length(binvec))) {
			temp = value
			value = temp*2+as.integer(binvec[i])
		}
		
	}
	return(as.integer(value))
}

xyzt2char <-
function(char) 
#convert a xyzt_unit raw element to space and time
{
	
	if(!is.character(char)) return(paste('xyzt_units not decipherable.'))
	
	bitvec = rawToBits(charToRaw(char))
	space = rev(bitvec[1:3])
	time = rev(bitvec[4:6])
	
	space_int = bin2dec(space)
	time_int = bin2dec(time)*8
	
	space_char = 'unknown'
	if(space_int==1) space_char = 'meters'
	if(space_int==2) space_char = 'millimeters'
	if(space_int==3) space_char = 'micrometers'
	
	int = 'unknown'
	if(time_int==8)  time_char = 'seconds'
	if(time_int==16) time_char = 'milliseconds'
	if(time_int==24) time_char = 'microseconds'
	if(time_int==32) time_char = 'Hertz'
	if(time_int==40) time_char = 'ppm'
	if(time_int==48) time_char = 'rad/sec'
	
	return(paste('xyz in ',space_char,' | time in ',time_char,'\n',sep=''))

}

removeNullString <- function(raw) 
#remove embedded and trailing nullstrings
{
	rm = which(raw==0)
	
	if(length(rm)>0) raw = raw[-rm]
	
	return(raw)
	
}

