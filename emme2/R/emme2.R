#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script to read and write to the EMME/2 databank
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Created: Ben Stabler 2/16/03 benjamin.stabler@odot.state.or.us
# Updated: Ben Stabler 6/11/03
# Updated: Ben Stabler 6/17/03
# Updated: Ben Stabler 8/3/04

# Updated 1/14/2013 to replace real() with double(), added NAMESPACE,
#  updated documentation since doesn't work for EMME/4 banks, and
#  added Peter Schmiedeskamp's functions 

# Copyright (C) 2002  Oregon Department of Transportation
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTION TO READ DATABANK FILE OFFSETS FOR READING THE DATABANK
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# File 0 consists of the offsets and number of records and bytes per file
read.file0 <- function(bank) {
	#bank is a string of the file name
	#Returns databank file offets (similar to EMME/2 module 1.1.5)
			
	#Open the databank for binary reading	
	infile <- file(bank, "rb")
	
	#Read in the file type (!=99 then EMME/2ban was created with EMME/2 version 9.x)
	file.type <- readBin(infile, integer(), 1, 4)
	
	if (file.type==99) {
		#OLD DATABANK STRUCTURE
		# The old databank structure uses bit packing and is
		# more complicated than the new databank structure.
		
		# Function to convert a long byte decimal value into bit form
		declong2bin <- function(decimal) {
			#Convert a long byte decimal value into bit form
			remainder <- NULL
			while (decimal>1) {
				remainder <- c(decimal%%2, remainder)
				decimal <- trunc(decimal/2)
			}
			if (decimal==1) remainder <- c(decimal, remainder) 
			pad <- 32-length(remainder)
			remainder <- c(rep(0, pad), remainder)
			remainder <- rev(remainder)
			remainder
		}

		# Function to convert a binary number to integer format
		bin2dec <- function(bin.vector) {
			iterations <- length(bin.vector)
			times <- 2^seq(0,31)
			number <- 0
			for (i in 1:iterations) {
				number <- number + bin.vector[i]*times[i]
			}
			number
		}
		
		# Read in all the words and convert to readable format	
		file0 <- NULL
		for (i in 1:99) {
			seek(infile, where=i*4, origin="start")
			word1 <- readBin(infile,integer(), 1, 4)
			seek(infile, where=(i+100)*4, origin="start")
			word2 <- readBin(infile,integer(), 1, 4)
			#Convert to binary format
			word1b<-declong2bin(word1)
			word2b<-declong2bin(word2)
			#Subset bits for databank properties (bit unpack)
			offsetb <- word1b[c(1:28,32)]
			typeb <- word1b[c(29,30)]
			#Subset bits for databank properties (bit unpack)
			wordrecb <- word2b[c(1:21,32,31)]
			recb <- word2b[c(22:31)]
			#Convert binary to integer
			offset <- bin2dec(offsetb)
			type <- bin2dec(typeb)
			wordrec <- bin2dec(wordrecb)
			rec <- bin2dec(recb)
			#Concatenate databank dimension properties and rbind to file0 
			file.data <- c(offset, rec, wordrec, type)
			file0 <- rbind(file0, file.data)
			}	
				
	} else {	
		#NEW DATABANK STRUCTURE
		file0.data <- readBin(infile, integer(), 400, 4)
		file0 <- cbind(file0.data[1:99],file0.data[101:199],file0.data[201:299],file0.data[301:399])
	}
	
	close(infile)
	colnames(file0) <- c("offset","records","words/rec","type")
	rownames(file0) <- 1:99	
	#Return databank file offets (similar to EMME/2 module 1.1.5)
	file0
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTION TO READ THE EMME/2 DATABANK FILE 1 INFORMATION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# File 1 consists of the global and scenario parameters
read.file1 <- function(bank, file0) {
	#bank is a string of the file name
	#file0 is the databank metadata data frame created with read.file0()
			
	#Open the databank for binary reading	
	infile <- file(bank, "rb")
	
	#Seek to global and scenario parameter file 1 position
	gsp.offset <- file0[1,1]
	seek(infile, where=gsp.offset*4, origin="start")
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# GLOBAL PARAMETERS
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	#Define return list to save global parameters
	return <- list()
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Read in the global parameters (file1)
	file1 <- readBin(infile, integer(), 80, 4)
	
	#mscen - Maximum number of scenarios
	#mcent - Maximum number of centroids
	#mnode - Maximum number of nodes
	#mlink - Maximum number of links
	#mturn - Maximum number of turn penalty tables
	#mline - Maximum number of transit lines
	#mlseg - Maximum total number of line segments
	#mmat  - Maximum number of matrices
	#mfunc - Maximum number of functions per class
	#moper - Maximum total number of operators for all functions class
	
	names(file1) <- c("ldi","ldo","lgi","lgo","ldai","ldao",
	"lero","llio","lrep","lgraph","iphys1","iphys2","iphys3",
	"iphys4","iphys5","iphys6","iphys7","iphys8","iphys9","iphys10",
	"kmod","idev","ishort","lpsiz","ipge","idat","iusr","itpter",
	"itppri","itpplo","nexdg","nlerr","igcmd","modsid","iscen",
	"imodl","lmodl","icgm","imfb","ierop","klu","kcu","keu","iscpu",
	"larrow","blank","blank","blank","blank","idbrel","mscen","mcent",
	"mnode","mlink","mturn","mline","mlseg","mmat","mfunc","moper",
	rep("blank",20))
	
	#Add the global parameters to the return list
	return <- c(return, list(global=file1))
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# SCENARIO PARAMETERS
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
	#Define scenario parameter names	
	scenario.data.names <- c("ncent","nnode","nlink","nturn",
	"nline","nlseg","nkturn","istats","itsimp","blank","mpmat1",
	"mgauto","mgtran","mgadd","mgadt","blank","blank","blank",
	"blank","blank","blank","mtimau","mtimtr","mboatr","mwaitr",
	"mauxtr","minvtr","mgnatr","mnbotr","mw1tr","mcadt","mautoc",
	"mpmat4","mpmat2","mpmat3","mcadd","mindfa","mwpqau","mfpqau",
	"mpmat5","mpmat6","litau","lgapau","lepsau","iterau","istopc",
	"ixlmax","iaddop","iaddlu1","iaddlu2","itsau","littr","modtra",
	"itimtr","iwtf","iwtw","iatw","ittw","lefhdw","modimp","itstr",
	"blank","npauto","nvauto","nvassc","nvdadc","nvadda","nvtrac",
	"blank","blank","blank","iadtop","iadtlu1","iadtlu2","iadtat1",
	"iadtat2","iadtat3","iadtat4","blank","blank")
	
	#Read in scenario data
	for (i in 1:return$global["mscen"]) {
		scenario.data <- readBin(infile, integer(), 80, 4)
		names(scenario.data) <- scenario.data.names
		return <- c(return, list(i=scenario.data))
	}
	
	#Close the file connection
	close(infile)
	
	#Return a named list of global parameters and scenario parameters
	names(return) <- c("global",letters[seq(1,return$global["mscen"],1)])
	return
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTION TO READ THE MATRIX DIRECTORY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read.matdir <- function(bank, file0, mmat) {
	#bank is a string of the file name
	#file0 is the databank metadata data frame created with read.file0()
	#mmat is the maximum number of matrices defined for the bank
	# and is created by the read.file1 function
			
	#Open the databank for binary reading to find byte stream seek values
	infile <- file(bank, "rb")

	#Seek to matrix directory file 60 position
	mat.offset <- file0[60,1]
	seek(infile, where=mat.offset*4, origin="start")
	
	# Matrix directory (file60)
	matrix.types <- c("ms","mo","md","mf")
	
	cflag <- list()
	for (type in matrix.types) {	
		cflag.temp <- readBin(infile, integer(), 4*mmat, 1)
		cflag.temp <- matrix(cflag.temp, , 4, byrow=T)
		colnames(cflag.temp) <- c("defined","columnwise","read-only","futruse")
		cflag <- c(cflag, list(cflag.temp))
	}
	names(cflag) <- matrix.types
	#Cflag is a 1 bit bit pattern of four values comprising one byte
	#Cflag [,1] value of 0 = matrix not defined, value of 1 = matrix defined
	
	mat.time <- list()
	for (type in matrix.types) {	
		mat.time.stamp <- readBin(infile, integer(), mmat, 4)
		mat.time <- c(mat.time, list(mat.time.stamp))
	}
	names(mat.time) <- matrix.types	
	#mat.time is a time stamp for each matrix
	
	mat.name <- list()
	for (type in matrix.types) {	
		mat.name.type <- NULL
		for (i in 1:mmat) {
			mat.name.temp <- readChar(infile, 12)
			mat.name.temp <- gsub(" +$","",mat.name.temp)
			mat.name.temp <- gsub("  ","",mat.name.temp)
			mat.name.type <- c(mat.name.type, mat.name.temp)
		}
		mat.name <- c(mat.name, list(mat.name.type))
	}
	names(mat.name) <- matrix.types	
	#mat.name is the matrix name
	
	mat.desc <- list()
	for (type in matrix.types) {	
		mat.desc.type <- NULL
		for (i in 1:mmat) {
			mat.desc.temp <- readChar(infile, 80)
			mat.desc.temp <- gsub(" +$","",mat.desc.temp)
			mat.desc.temp <- gsub("  ","",mat.desc.temp)
			mat.desc.type <- c(mat.desc.type, mat.desc.temp)
		}
		mat.desc <- c(mat.desc, list(mat.desc.type))
	}
	names(mat.desc) <- matrix.types	
	#mat.desc is the matrix description
	
	#Create mat.dir results of matrix directory
	mat.dir <- list()
	mat.dir <- c(mat.dir, list(ms=cbind(name=mat.name$ms, desc=mat.desc$ms)))
	mat.dir <- c(mat.dir, list(mo=cbind(name=mat.name$mo, desc=mat.desc$mo)))
	mat.dir <- c(mat.dir, list(md=cbind(name=mat.name$md, desc=mat.desc$md)))
	mat.dir <- c(mat.dir, list(mf=cbind(name=mat.name$mf, desc=mat.desc$mf)))

	#Close the file connection and return the matrix directory
	close(infile)
	mat.dir
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTION TO READ ALL MSs FROM THE DATABANK
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read.ms <- function(bank, file0) {
	#bank is a string of the file name
	#file0 is the databank metadata data frame created with read.file0()
	
	#Open the databank for binary reading to find byte stream seek values
	infile <- file(bank, "rb")

	#Seek to origin matrix file 61 position
	mat.offset <- file0[61,1]
	seek(infile, where=mat.offset*4, origin="start")
	
	#Read in the matrix data
	ms <- readBin(infile, double(), file0[61,3], 4)
	names(ms) <- 1:length(ms)
		
	#Close the infile connection and return the matrix
	close(infile)
	ms
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTION TO READ A MO FROM THE DATABANK
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read.mo <- function(numname, bank, file0, mcent, mat.dir) {
	#bank is a string of the file name
	#numname is the mf number or name as a string to read in
	#file0 is the databank metadata data frame created with read.file0()
	#mcent is the maximum number of centroids defined for the bank
	# and is created by the read.file1 function
	#mat.dir is the matrix directory object created by read.matdir()
	
	#Open the databank for binary reading to find byte stream seek values
	infile <- file(bank, "rb")

	#Seek to origin matrix file 62 position
	mat.offset <- file0[62,1]
	seek(infile, where=mat.offset*4, origin="start")
	
	#Lookup matrix number from name
	if (is.character(numname)) {
		numname <- gsub(" +$","",numname)
		number <- which(mat.dir$mo[,1]==numname)
		if (length(number)==0) { stop("Matrix Not Found") }
	} else { number <- numname }
	
	#Seek to the specific matrix in the origin matrix file
	specific.offset <- (number-1)*mcent
	seek(infile, where=specific.offset*4, origin="current")
	
	#Read in the matrix data
	mo <- readBin(infile, double(), mcent, 4)
		
	#Close the infile connection and return the matrix
	close(infile)
	mo
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTION TO READ A MD FROM THE DATABANK
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read.md <- function(numname, bank, file0, mcent, mat.dir) {
	#bank is a string of the file name
	#numname is the mf number or name as a string to read in
	#file0 is the databank metadata data frame created with read.file0()
	#mcent is the maximum number of centroids defined for the bank
	# and is created by the read.file1 function
	#mat.dir is the matrix directory object created by read.matdir()
	
	#Open the databank for binary reading to find byte stream seek values
	infile <- file(bank, "rb")

	#Seek to destination matrix file 63 position
	mat.offset <- file0[63,1]
	seek(infile, where=mat.offset*4, origin="start")
	
	#Lookup matrix number from name
	if (is.character(numname)) {
		numname <- gsub(" +$","",numname)
		number <- which(mat.dir$md[,1]==numname)
		if (length(number)==0) { stop("Matrix Not Found") }
	} else { number <- numname }
		
	#Seek to the specific matrix in the destination matrix file
	specific.offset <- (number-1)*mcent
	seek(infile, where=specific.offset*4, origin="current")
	
	#Read in the matrix data
	md <- readBin(infile, double(), mcent, 4)
		
	#Close the infile connection and return the matrix
	close(infile)
	md
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTION TO READ A MF FROM THE DATABANK
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read.mf <- function(numname, bank, file0, mcent, mat.dir) {
	#bank is a string of the file name
	#numname is the mf number or name as a string to read in
	#file0 is the databank metadata data frame created with read.file0()
	#mcent is the maximum number of centroids defined for the bank
	# and is created by the read.file1 function
	#mat.dir is the matrix directory object created by read.matdir()
	
	#Open the databank for binary reading to find byte stream seek values
	infile <- file(bank, "rb")

	#Seek to full matrix file 64 position
	mat.offset <- file0[64,1]
	seek(infile, where=mat.offset*4, origin="start")
	
	#Lookup matrix number from name
	if (is.character(numname)) {
		numname <- gsub(" +$","",numname)
		number <- which(mat.dir$mf[,1]==numname)
		if(length(number)>1) {print(paste("Warning, non-unique matrix. Using the first: ", numname))}
		number <- number[1]
		if (length(number)==0) { stop("Matrix Not Found") }
	} else { number <- numname }
	
	#Seek to the specific matrix in the full matrix file
	#+++++++++++++++++++++++++++++++++
	#Start Steve Hansen Edit
	#+++++++++++++++++++++++++++++++++
	if(number>1){
		mf.offset <- mcent*mcent
		iterations <- number-1
		for (i in 1:iterations){
			seek(infile, where=mf.offset*4, origin="current")
		}
	}
	#+++++++++++++++++++++++++++++++++
	#End Steve Hansen Edit
	#+++++++++++++++++++++++++++++++++
	#Read in the matrix data
	mf.temp <- readBin(infile, double(), mcent*mcent, 4)
	mf <- matrix(mf.temp, mcent, mcent, byrow=T)
		
	#Close the infile connection and return the matrix
	close(infile)
	mf
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTION TO WRITE A MF TO THE DATABANK
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.mf <- function(data, numname, bank, file0, mcent, mmat, mat.dir, newname=NULL, newdesc=NULL) {
	#data is either a vector or matrix
	#numname is the mf number or name of the existing matrix to replace
	#bank is a string of the file name
	#file0 is the databank metadata data frame created with read.file0()
	#Note that EMME/2 stores mfs in row-major order
	#mcent is the maximum number of centroids defined for the bank
	# and is created by the read.file1 function
	#mmat is the maximum number of matrices allowed for each type
	#mat.dir is the matrix directory object created by read.matdir()
	#newname is the new name of the matrix to write out
	#newdesc is the new description of the matrix to write out
				
	#Open the databank for binary reading to find byte stream seek values
	outfile <- file(bank, "r+b")

	#Seek to full matrix file 64 position
	mat.offset <- file0[64,1]
	seek(outfile, where=mat.offset*4, origin="start")
	
	#Lookup matrix number from name
	if (is.character(numname)) {
		numname <- gsub(" +$","",numname)
		number <- which(mat.dir$mf[,1]==numname)
		if (length(number)==0) { stop("Matrix Not Found") }
	} else { number <- numname }
			
	#Seek to the specific matrix in the full matrix file
	#+++++++++++++++++++++++++++++++++
	#Start Steve Hansen Edit
	#+++++++++++++++++++++++++++++++++
	if(number>1){
		mf.offset <- mcent*mcent
		iterations <- number-1
		for (i in 1:iterations){
			seek(outfile, where=mf.offset*4, origin="current", rw="write")
		}
	}
	#+++++++++++++++++++++++++++++++++
	#End Steve Hansen Edit

	#If data is a matrix, then convert the matrix to vector form
	if (is.matrix(data)) {
		data <- as.vector(t(data))
	}
	
	#Write mfnumber to the databank
	writeBin(as.double(data), outfile, 4)
		
	#Seek to matrix directory file 60 position
	mat.offset <- file0[60,1]
	seek(outfile, where=mat.offset*4, origin="start", rw="write")
	#Seek to mf part of cflag
	seek(outfile, where=mmat*4*3, origin="current", rw="write")
	#Seek to cflag entry for specifc mf and write 1 to tag matrix as defined
	seek(outfile, where=(number-1)*4, origin="current", rw="write")
	writeBin(as.integer(1), outfile, 1)
	
	if (!is.null(newname)) {
		#Seek to matrix directory file 60 position
		seek(outfile, where=mat.offset*4, origin="start", rw="write")
		#Seek to name part of matrix directory
		seek(outfile, where=(mmat*4*4)+(mmat*4*4), origin="current", rw="write")
		
		#Seek to name entry for mfs
		seek(outfile, where=mmat*4*3*3, origin="current", rw="write")
		#Seek to specific name entry for matrix
		seek(outfile, where=(number-1)*4*3, origin="current", rw="write")
		
		#Create name format (2 chars 2 spaces 2 chars 2 spaces up to 12)
		newname <- substring(newname,c(1,3,5),c(2,4,6))
		newname <- paste(newname, collapse="  ")
		newname <- paste(newname, paste(rep(" ",12-nchar(newname)), collapse=""), collapse="")
		#Write matrix name (no spaces allowed)
		writeChar(as.character(newname), outfile, 12, eos=NULL)
	}
	
	if (!is.null(newdesc)) {
		#Seek to matrix directory file 60 position
		seek(outfile, where=mat.offset*4, origin="start", rw="write")
		#Seek to description part of matrix directory
		seek(outfile, where=(mmat*4*4)+(mmat*4*4)+(mmat*4*4*3), origin="current", rw="write")
		
		#Seek to description entry for mfs
		seek(outfile, where=mmat*4*20*3, origin="current", rw="write")
		#Seek to specific description entry for matrix
		seek(outfile, where=(number-1)*4*20, origin="current", rw="write")
		
		#Create description format (2 chars 2 spaces 2 chars 2 spaces up to 80)
		newdesc <- substring(newdesc,seq(1,40,2),seq(2,40,2))
		newdesc <- paste(newdesc, collapse="  ")
		newdesc <- paste(newdesc, paste(rep(" ",80-nchar(newdesc)), collapse=""), collapse="")
		#Write matrix description (spaces allowed)
		writeChar(as.character(newdesc), outfile, 80, eos=NULL)
	}
	
	#+++++++++++++++++++++++++++++++++
	#Start Steve Hansen Edit
	#+++++++++++++++++++++++++++++++++
	
	## Write datestamp ##
	#Seek to matrix directory file 60 position
	mat.offset <- file0[60,1]
	seek(outfile, where=mat.offset*4, origin="start", rw="write")
	#Seek to part of matrix directory file that stores the timestamp for matrix_number
	seek(outfile, where=mmat*28+(numname-1)*4, origin="current", rw="write")
	emme2time <- get.emme2.time(Sys.time())
	writeBin(as.integer(emme2time), outfile, 4)	
	#+++++++++++++++++++++++++++++++++
	#End Steve Hansen Edit
	#+++++++++++++++++++++++++++++++++
	
	#Close the file connection
	close(outfile)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTION TO READ LINK SPEED, CAPACITY, AND VDF FROM A SCENARIO
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# First need to call read.file0 and read.file1 
# Brian Gregor, Brian.J.GREGOR@odot.state.or.us

read.link.data <- function (bank, scen.num, file0, mscen, mlink, mnode){

    # Identify the file to read
    infile <- file(bank, "rb")
    
    # Read the node data
    seek(infile, where = file0[6, 1] * 4, origin = "start")
    seek(infile, where = (scen.num - 1) * mnode * 4, origin = "current")
    node.data <- readBin(infile, integer(), mnode, 4)

    # Read the pointer to j node
    seek(infile, where = file0[9, 1] * 4, origin = "start")
    seek(infile, where = (scen.num - 1) * mnode * 4, origin = "current")
    pointer.to.j.node <- readBin(infile, integer(), mnode, 4)
    pointer.to.j.node <- diff(pointer.to.j.node)

    # Read the j node
    seek(infile, where = file0[11, 1] * 4, origin = "start")
    seek(infile, where = (scen.num - 1) * mlink * 4, origin = "current")
    j.node <- readBin(infile, integer(), mlink, 4)

    # Read the link length
    seek(infile, where = file0[12, 1] * 4, origin = "start")
    seek(infile, where = (scen.num - 1) * mlink * 4, origin = "current")
    link.length <- readBin(infile, integer(), mlink, 4)
    link.length <- link.length/100

    # Read the link type
    seek(infile, where = file0[14, 1] * 4, origin = "start")
    seek(infile, where = (scen.num - 1) * mlink * 4, origin = "current")
    link.type <- readBin(infile, integer(), mlink, 4)

    # Read the vdf and number of lanes
    seek(infile, where = file0[15, 1] * 4, origin = "start")
    seek(infile, where = (scen.num - 1) * mlink * 4, origin = "current")
    link.lanes.vdf <- readBin(infile, integer(), mlink, 4)
    link.vdf <- as.numeric(substring(link.lanes.vdf, 1, 1))
    link.lanes <- as.numeric(substring(link.lanes.vdf, 2, 3))/10

    # Read ul1 and ul2
    seek(infile, where = file0[16, 1] * 4, origin = "start")
    seek(infile, where = (scen.num - 1) * mlink * 4, origin = "current")
    link.ul1 <- readBin(infile, double(), mlink, 4)
    seek(infile, where = (file0[16, 1] * 4 + mscen * mlink * 
        4), origin = "start")
    seek(infile, where = (scen.num - 1) * mlink * 4, origin = "current")
    link.ul2 <- readBin(infile, double(), mlink, 4)

    # Read timau
    seek(infile, where = file0[17, 1] * 4, origin = "start")
    seek(infile, where = (scen.num - 1) * mlink * 4, origin = "current")
    timau <- readBin(infile, double(), mlink, 4)
    
    # Read volau
    seek(infile, where = file0[18, 1] * 4, origin = "start")
    seek(infile, where = (scen.num - 1) * mlink * 4, origin = "current")
    volau <- readBin(infile, double(), mlink, 4)
    
    close(infile)
    list(node.data, pointer.to.j.node, j.node, length = link.length, 
        type = link.type, vdf = link.vdf, lanes = link.lanes, 
        ul1 = link.ul1, ul2 = link.ul2, timau = timau, volau = volau)
    }



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTION TO READ NODE DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# First need to call read.file0 and read.file1 
read.nodes <- function(bank, scen.num, file0, mscen, mlink, mnode) {
	#bank is a string of the file name
	#scen.num is the scenario number to read from (in EMME/2 order - not named number)
	#file0 is the databank metadata data frame created with read.file0()
	#mscen is the maximum number of scenarios defined for the bank
	#mlink is the maximum number of links defined for the bank
	#mnode is the maximum number of nodes defined for the bank
			
	infile<-file(bank, "rb")
	
	#Read in node data
	seek(infile, where=file0[6,1]*4, origin="start")
	seek(infile, where=(scen.num-1)*mnode*4, origin="current")
	node.data <- readBin(infile,integer(),mnode,4)
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Read in X and Y coordinates of nodes (file7)
	seek(infile, where=file0[7,1]*4, origin="start")
	seek(infile, where=(scen.num-1)*mnode*4, origin="current")
	x <- readBin(infile,double(), mnode, 4)
	
	seek(infile, where=file0[7,1]*4, origin="start")
	seek(infile, where=(scen.num-1)*mnode*4+mscen*mnode*4, origin="current")
	y <- readBin(infile,double(), mnode, 4)
	
	close(infile)
	nodes <- data.frame(id=node.data, x=x, y=y)
	nodes <- nodes[nodes$id>0,]
	nodes
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTION TO PLOT THE BASE NETWORK
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plotLinks <- function (tofrom, nodes, title, ...) 
{
    plot(nodes[, 2], nodes[, 3], type = "n", xlab = "X", ylab = "Y", 
        main = title)
    fnode.xy <- nodes[match(tofrom[, 1], nodes[, 1]), ]
    tnode.xy <- nodes[match(tofrom[, 2], nodes[, 1]), ]
    ftxy <- cbind(fnode.xy, tnode.xy)
    invisible(apply(ftxy, 1, function(x) lines(x[c(2, 5)], x[c(3, 
        6)], pch = 20, type = "o", ...)))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTION TO BUILD A FNODE TNODE TABLE 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To build the I J links
ftnode <- function(node.data, outgoing.links, jnode, mlink) {
	outgoing.links[outgoing.links<0] <- 0
	outgoing.links <- c(outgoing.links,0)
	i.nodes <- rep(node.data,outgoing.links)
	i.nodes <- c(i.nodes, rep(0, mlink-length(i.nodes)))
	j.nodes <- node.data[jnode]
	j.nodes <- c(j.nodes, rep(0, mlink-length(j.nodes)))
	ijnode <- cbind(fnode=i.nodes,tnode=j.nodes)
	ijnode
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTION TO FORMAT A MATRIX FOR WRITING TO THE BANK 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Function to append zeros to matrix to write into databank
formatMf <- function(data, file1) {
     #Matrix build format is  m1 m2
     #                        m3 m4
     cols <- file1$global["mcent"] - ncol(data)
     rows <- file1$global["mcent"] - nrow(data)
     m1 <- data
     
     m2 <- matrix(0, nrow(m1), cols)
     m3 <- matrix(0, rows, ncol(m1))
     m4 <- matrix(0, rows, cols)
     rbind(cbind(m1, m2), cbind(m3,m4)) 
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTION TO CREATE AN INTEGER BASED ON THE EMME2 FORMULA FOR A DATESTAMP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get.emme2.time <- function(timestamp){
  #first split date from time
  date<-unlist(strsplit(as.character(timestamp)," "))[1]
  time<-unlist(strsplit(as.character(timestamp)," "))[2]
  year<-as.integer(unlist(strsplit(as.character(date),"-"))[1])
  month<-as.integer(unlist(strsplit(as.character(date),"-"))[2])
  day<-as.integer(unlist(strsplit(as.character(date),"-"))[3])
  hour<-as.integer(unlist(strsplit(as.character(time),":"))[1])
  minute<-as.integer(unlist(strsplit(as.character(time),":"))[2])
  second<-as.integer(unlist(strsplit(as.character(time),":"))[3])
  emme2time<-second+60*(minute+60*(hour+24*(day-1+31*(month-1+12*(year-1990)))))
  emme2time
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CONVENIENCE FUNCTIONS TO RETURN MATRICES AND DIRECTORIES IN DATA.FRAME FORMAT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Batch read a list of emme matrix names (short version), and return a merged data.frame
# DANGER: This function is not exactly memory-efficient. Don't attempt on low-memory systems without refactoring.
MFBatchFetch <- function(bank, matrixlist, useshortnames=FALSE) {
  # Fetch the first matrix off the list
  # (assumes first matrix is representative of the those following)
  print(matrixlist[1])
  
  if(useshortnames==TRUE){
    df <- MFFetch(bank, matrixlist[1],varlongname=matrixlist[1])
  }else{
    df <- MFFetch(bank, matrixlist[1])
  }
  
  # merging as opposed to cbind is slow, but enforces some sort of sanity checking
  for(m in matrixlist[-1]) {
    print(m)
    
    if(useshortnames==TRUE){
      new <- MFFetch(bank, m, varlongname=m)
    } else {
      new <- MFFetch(bank, m)
    }
    
    df <- merge(df, new, by=intersect(names(df), names(new)))
  }
  return(df)
}

# Return the named matrix as a neatly-formatted dataframe
MFFetch <- function(bank, matrixname, varlongname=NULL, valsonly=NULL) {
 
  library(reshape)

  file0 <- read.file0(bank)
  file1 <- read.file1(bank, file0)
  mat.dir <- read.matdir(bank, file0, file1$global["mmat"])
  mf <- read.mf(matrixname, bank, file0, file1$global["mcent"], mat.dir)
  mf <- melt(mf,)
  
  if(is.null(varlongname)) {
    # Convert this to a dataframe to avoid breaking my brain
    dirdf <- data.frame(mat.dir$mf, stringsAsFactors=FALSE)
    # Sometimes this won't be unique, hense the "[1]"
    varlongname <- dirdf[dirdf$name==matrixname, 2][1]
    varlongname <- gsub("\\s+", "_", varlongname)
    varlongname <- gsub("-", "_", varlongname)
    varlongname <- gsub("_+", "_", varlongname)
  }
  
  names(mf) <- c("orig", "dest", varlongname)
  
  if(is.null(valsonly)) {
    mf
  } else {
    mf[-c(1,2)]
  }
}

# Return a data.frame containing an emmebank's directory listing
MFDir <- function(bank) {
  # Boilerplate file index reading
  file0 <- read.file0(bank)
  file1 <- read.file1(bank, file0)
  # Return a more useful data type and filter out null entries
  mat.dir <- data.frame(read.matdir(bank, file0, file1$global["mmat"])$mf, stringsAsFactors=FALSE)
  mat.dir <- mat.dir[mat.dir$name != "", ]
  
  return(mat.dir)
}
