#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Shapefile Format - Read/Write shapefile format within R
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Created 2/13/03 Ben Stabler benjamin.stabler@odot.state.or.us
# Revised 6/11/03 Ben Stabler benjamin.stabler@odot.state.or.us
# Revised 7/7/03  Ben Stabler benjamin.stabler@odot.state.or.us
# Revised 7/23/03 Ben Stabler benjamin.stabler@odot.state.or.us

# Revised 8/10/05 Ben Stabler, benstabler@yahoo.com
# Reading of shape type 13 and 15 from Don MacQueen, macq@llnl.gov
# Conversion of simple polygons to shapefile format 
#  	from Manuel Chirouze, Manuel.Chirouze@benfieldgroup.com
# Small bug fix to write.dbf to deal with 1 column white space correctly
# Revised 8/15/05 Ben Stabler, to use write.dbf and read.dbf from foreign library
#  Added dp function as well

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
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DESCRIPTION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#source("shapefiles.R")
#myshape <- read.shapefile("myshapewithnoextension")
#myshape$shp - the shp file 
#myshape$shx - the index file
#myshape$dbf - the dbf file

#Each of the three lists of the myshape contain a header and some data. So there is...
#myshape$shp$header (the header info) and myshape$shp$shp (the geographic data)
#myshape$shx$header and myshape$shx$index (the index data)
#myshape$dbf$header and myshape$dbf$dbf (the dbf stored as a R data.frame)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Read in a point, line or polygon SHP file
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read.shp <- function(shp.name) {
	
	infile<-file(shp.name,"rb")

	#Header
	file.code <- readBin(infile,integer(),1, endian="big")
	unused <- readBin(infile,integer(),5, endian="big")
	file.length <- readBin(infile,integer(),1, endian="big")
	file.version <- readBin(infile,integer(),1, endian="little")
	shape.type <- readBin(infile,integer(),1, endian="little")
	xmin <- readBin(infile,double(),1, endian="little")
	ymin <- readBin(infile,double(),1, endian="little")
	xmax <- readBin(infile,double(),1, endian="little")
	ymax <- readBin(infile,double(),1, endian="little")
	zmin <- readBin(infile,double(),1, endian="little")
	zmax <- readBin(infile,double(),1, endian="little")
	mmin <- readBin(infile,double(),1, endian="little")
	mmax <- readBin(infile,double(),1, endian="little")
	header <- list(file.code,file.length,file.version,shape.type,xmin,ymin,xmax,ymax,zmin,zmax,mmin,mmax)
	names(header) <- c("file.code","file.length","file.version","shape.type","xmin","ymin","xmax","ymax","zmin","zmax","mmin","mmax")
	rm(file.code,file.length,file.version,shape.type,xmin,ymin,xmax,ymax,zmin,zmax,mmin,mmax,unused)
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	#Point Shapefile
	if (header$shape.type==1) {
		record <- integer()
		shape.type <- integer()
		x <- double()
		y <- double()
		record.num <- 1
		while (length(record.num)!=0) {
			#Record Number	
			record.num <- readBin(infile,integer(),1, 4, endian="big")
			if (length(record.num)==0) break
			record <- c(record,record.num)
			#Content Length
			content.length <- readBin(infile,integer(),1, 4, endian="big")
			#Shape Type (should be 1)
			temp.type <- readBin(infile,integer(),1, 4, endian="little")
			shape.type <- c(shape.type, temp.type)
			if (temp.type==1) {
				#X and Y values
				x <- c(x,readBin(infile,double(),1, 8, endian="little"))
				y <- c(y,readBin(infile,double(),1, 8, endian="little"))
				#To read in extra data created by GeoMedia shapefile export
				if (content.length>10) {
					temp <- readBin(infile,integer(), ((content.length-10)/2), 4, endian="big")
				}
			}
			if (temp.type==0) {
				#Null shape type has no data
				x <- c(x,NA)
				y <- c(y,NA)			
			}	
		}	
		point.data <- cbind(record=record, x=x, y=y, shape.type=shape.type)
		close(infile)
	}
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	#PolyLine, Polygon, PolyLineZ, and PolygonZ Shapefile
	if (header$shape.type==3 || header$shape.type==5 || header$shape.type == 13 || header$shape.type == 15) {
		shape <- list()
		record.num <- 1
		while (length(record.num)!=0) {
			#Record Number	
			record.num <- readBin(infile,integer(),1, endian="big")
			if (length(record.num)==0) break
			#Content Length
			content.length <- readBin(infile,integer(),1, endian="big")	
			#Shape Type
			shape.type <- readBin(infile, integer(), 1, endian="little")
			if (shape.type==3 || shape.type==5) {
				box <- readBin(infile, double(), 4, endian="little")
				names(box) <- c("xmin","ymin","xmax","ymax")
				num.parts <- readBin(infile, integer(), 1, endian="little")
				num.points <- readBin(infile, integer(), 1, endian="little")
				parts <- readBin(infile, integer(), num.parts, endian="little")
				points <- readBin(infile, double(), num.points*2, endian="little")
				X <- points[seq(1,(num.points*2),by=2)]
				Y <- points[seq(2,(num.points*2),by=2)]
				points <- data.frame(X,Y)
			}
			
			#Type 13 polylineZ (and 15 polygonZ since it is the same as 13)
			if (shape.type==13 || shape.type==15) {
			        box <- readBin(infile, double(), 4, endian = "little")
			        names(box) <- c("xmin", "ymin", "xmax", "ymax")
			        num.parts <- readBin(infile, integer(), 1, endian = "little")
			        num.points <- readBin(infile, integer(), 1, endian = "little")
			        parts <- readBin(infile, integer(), num.parts, endian = "little")
			        points <- readBin(infile, double(), num.points * 2, endian = "little")
			        X <- points[seq(1, (num.points * 2), by = 2)]
			        Y <- points[seq(2, (num.points * 2), by = 2)]
			        #Additional information for type 13 and 15
			        zrange <- readBin(infile,double(),2,endian='little')
			        Z <- readBin(infile,double(),num.points,endian='little')
			        mrange <- readBin(infile,double(),2,endian='little')
			        M <- readBin(infile,double(),num.points,endian='little')
			        points <- data.frame(X, Y, Z, M)
			}
					
			if (shape.type==0) {
				#Null shape type has no data
				box <- rep(NA,4)
				num.parts <- NA
				num.points <- NA
				parts <- NA
				points <- data.frame(X=NA,Y=NA)				
			}	
			shape.info <- list(record=record.num, content.length=content.length, shape.type=shape.type, box=box, num.parts=num.parts, num.points=num.points, parts=parts, points=points)
			shape <- c(shape,list(shape.info))
		}
			close(infile)
	}
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	if (header$shape.type==1) { 
		return.val <- point.data 
	}
	if (header$shape.type==3 || header$shape.type==5 || header$shape.type == 13 || header$shape.type == 15) { 
		return.val <- shape
	}
	list(shp=return.val,header=header)
}	


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Write out a point, line or polygon SHP file
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write.shp <- function(shp, out.name) {
	
	outfile<-file(out.name,"wb")
	
	#Header
	writeBin(as.integer(9994), outfile, endian="big")
	writeBin(as.integer(rep(0,5)), outfile, endian="big")
	writeBin(as.integer(shp$header$file.length), outfile, endian="big")
	writeBin(as.integer(1000), outfile, endian="little")
	writeBin(as.integer(shp$header$shape.type), outfile, endian="little")
	writeBin(as.double(shp$header$xmin), outfile, endian="little")
	writeBin(as.double(shp$header$ymin), outfile, endian="little")
	writeBin(as.double(shp$header$xmax), outfile, endian="little")
	writeBin(as.double(shp$header$ymax), outfile, endian="little")
	writeBin(as.double(shp$header$zmin), outfile, endian="little")
	writeBin(as.double(shp$header$zmax), outfile, endian="little")
	writeBin(as.double(shp$header$mmin), outfile, endian="little")
	writeBin(as.double(shp$header$mmax), outfile, endian="little")
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	if (shp$header$shape.type==1) {
		for (record in 1:nrow(shp$shp)) {
			if (shp$shp[record,4]==1) {
				#Record Number
				writeBin(as.integer(shp$shp[record,1]), outfile, endian="big")
				#Content Length
				writeBin(as.integer(10), outfile, endian="big")
				#Shape Type
				writeBin(as.integer(1), outfile, endian="little")
				#X - For Points only
				writeBin(as.double(shp$shp[record,2]), outfile, endian="little")
				#Y - For Points only
				writeBin(as.double(shp$shp[record,3]), outfile, endian="little")
			}
			if (shp$shp[record,4]==0) {
				#Record Number
				writeBin(as.integer(shp$shp[record,1]), outfile, endian="big")
				#Content Length
				writeBin(as.integer(4), outfile, endian="big")
				#Shape Type
				writeBin(as.integer(0), outfile, endian="little")
			}
		}
	}
	
	if (shp$header$shape.type==3 || shp$header$shape.type==5 || shp$header$shape.type == 13 || shp$header$shape.type == 15) {
		for (record in 1:length(shp$shp)) {
			if (shp$shp[[record]]$shape.type==3 || shp$shp[[record]]$shape.type==5) {
				#Record Number
				writeBin(as.integer(shp$shp[[record]]$record), outfile, endian="big")
				#Content Length
				writeBin(as.integer(shp$shp[[record]]$content.length), outfile, endian="big")
				#Shape Type
				writeBin(as.integer(shp$shp[[record]]$shape.type), outfile, endian="little")
				#Bounding Box
				writeBin(as.double(shp$shp[[record]]$box[1]), outfile, endian="little")
				writeBin(as.double(shp$shp[[record]]$box[2]), outfile, endian="little")
				writeBin(as.double(shp$shp[[record]]$box[3]), outfile, endian="little")
				writeBin(as.double(shp$shp[[record]]$box[4]), outfile, endian="little")
				#Number of parts (segments)
				writeBin(as.integer(shp$shp[[record]]$num.parts), outfile, endian="little")
				#Number of total points
				writeBin(as.integer(shp$shp[[record]]$num.points), outfile, endian="little")
				#Parts
				writeBin(as.integer(shp$shp[[record]]$parts), outfile, endian="little")
				#Points - merge the X and Y points into one vector
				point.stream <- rep(0,length(shp$shp[[record]]$points$X)*2)
				point.stream[seq(1,length(shp$shp[[record]]$points$X)*2,by=2)] <- shp$shp[[record]]$points$X
				point.stream[seq(2,length(shp$shp[[record]]$points$Y)*2,by=2)] <- shp$shp[[record]]$points$Y
				writeBin(as.double(point.stream), outfile, endian="little")
				#Need to check to make sure first and last vertex are the same for Polygon
			}
			
			#Type 13 polylineZ (and 15 polygonZ since it is the same as 13)
			if (shp$shp[[record]]$shape.type==13 || shp$shp[[record]]$shape.type==15) {
			
				#Record Number
				writeBin(as.integer(shp$shp[[record]]$record), outfile, endian="big")
				#Content Length
				writeBin(as.integer(shp$shp[[record]]$content.length), outfile, endian="big")
				#Shape Type
				writeBin(as.integer(shp$shp[[record]]$shape.type), outfile, endian="little")
				#Bounding Box
				writeBin(as.double(shp$shp[[record]]$box[1]), outfile, endian="little")
				writeBin(as.double(shp$shp[[record]]$box[2]), outfile, endian="little")
				writeBin(as.double(shp$shp[[record]]$box[3]), outfile, endian="little")
				writeBin(as.double(shp$shp[[record]]$box[4]), outfile, endian="little")
				#Number of parts (segments)
				writeBin(as.integer(shp$shp[[record]]$num.parts), outfile, endian="little")
				#Number of total points
				writeBin(as.integer(shp$shp[[record]]$num.points), outfile, endian="little")
				#Parts
				writeBin(as.integer(shp$shp[[record]]$parts), outfile, endian="little")
				
				#Points - merge the X and Y points into one vector
				point.stream <- rep(0,length(shp$shp[[record]]$points$X)*2)
				point.stream[seq(1,length(shp$shp[[record]]$points$X)*2,by=2)] <- shp$shp[[record]]$points$X
				point.stream[seq(2,length(shp$shp[[record]]$points$Y)*2,by=2)] <- shp$shp[[record]]$points$Y
				writeBin(as.double(point.stream), outfile, endian="little")
				#Need to check to make sure first and last vertex are the same for Polygon
			
				#Write out Z Range and Z values
				minZ <- min(shp$shp[[record]]$points$Z)
				maxZ <- max(shp$shp[[record]]$points$Z)
				writeBin(as.double(minZ), outfile, endian='little')
				writeBin(as.double(maxZ), outfile, endian='little')
				writeBin(as.double(shp$shp[[record]]$points$Z), outfile, endian='little')
				
				#Write out M Range and M values
				minM <- min(shp$shp[[record]]$points$M)
				maxM <- max(shp$shp[[record]]$points$M)
				writeBin(as.double(minM), outfile, endian='little')
				writeBin(as.double(maxM), outfile, endian='little')				
				writeBin(as.double(shp$shp[[record]]$points$M), outfile, endian='little')
			}
			
			
			if (shp$shp[[record]]$shape.type==0) {
				#Record Number
				writeBin(as.integer(shp$shp[[record]]$record), outfile, endian="big")
				#Content Length
				writeBin(as.integer(shp$shp[[record]]$content.length), outfile, endian="big")
				#Shape Type
				writeBin(as.integer(shp$shp[[record]]$shape.type), outfile, endian="little")
			}
		}
	}
	close(outfile)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Read in the SHX file (index)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

read.shx <- function(shx.name) {
	
	infile<-file(shx.name,"rb")
	
	#Header
	file.code <- readBin(infile,integer(),1, endian="big")
	unused <- readBin(infile,integer(),5, endian="big")
	file.length <- readBin(infile,integer(),1, endian="big")
	file.version <- readBin(infile,integer(),1, endian="little")
	shape.type <- readBin(infile,integer(),1, endian="little")
	xmin <- readBin(infile,double(),1, endian="little")
	ymin <- readBin(infile,double(),1, endian="little")
	xmax <- readBin(infile,double(),1, endian="little")
	ymax <- readBin(infile,double(),1, endian="little")
	zmin <- readBin(infile,double(),1, endian="little")
	zmax <- readBin(infile,double(),1, endian="little")
	mmin <- readBin(infile,double(),1, endian="little")
	mmax <- readBin(infile,double(),1, endian="little")
	header <- list(file.code,file.length,file.version,shape.type,xmin,ymin,xmax,ymax,zmin,zmax,mmin,mmax)
	names(header) <- c("file.code","file.length","file.version","shape.type","xmin","ymin","xmax","ymax","zmin","zmax","mmin","mmax")
	records <- (file.length-50)/4
	
	#Record Data
	index.data <- matrix(0,records,2)
	colnames(index.data) <- c("Offset","Length")
	for (record in 1:records) {
		#First record offset is 50
		index.data[record,1] <- readBin(infile,integer(),1, endian="big")
		#Content length
		index.data[record,2] <- readBin(infile,integer(),1, endian="big")
	}
	close(infile)
	list(index=index.data,header=header)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to calculate header information for a point file
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Used to clean up GeoMedia shapefile export files

calc.header <- function(shapefile) {
	#1) Calculate shp header info
	#Calculate number of bytes for null record headers
	num.null.record.bytes <- length(shapefile$shp$shp[,2][is.na(shapefile$shp$shp[,2])]) * 4 
	#Calculate number of bytes for complete record record headers
	num.record.header.bytes <- length(shapefile$shp$shp[,2]) * 8
	#Calculate number of bytes for records (both geoMedia and generic shapefile)
	#Should be value 20 (56 for GeoMedia files)
	num.record.bytes <- (length(shapefile$shp$shp[,2]) - length(shapefile$shp$shp[,2][is.na(shapefile$shp$shp[,2])])) * 20
	#Sum the bytes to get total bytes (100 = header bytes)
	total.bytes <- 100 + sum(num.null.record.bytes, num.record.header.bytes, num.record.bytes)
	#Divide by 2 to get number of 16-bit words 
	file.length <- total.bytes / 2
	#Replace file.length in shapefile header with new file.length calculation
	shapefile$shp$header$file.length <- file.length
	#Replace lengths of 28 with 10 for shx from GeoMedia
	
	#2) Calculate shx info
	if (shapefile$shx$index[1,2]==28) {
		shapefile$shx$index[,2][shapefile$shx$index[,2]==28] <- 10
	}
	#Set the content.length to 14 in order to calculate offset
	shapefile$shx$index[,2][shapefile$shx$index[,2]==10]<-14
	#Set null value content.lengths from 2 to 6 to calculate offsets
	shapefile$shx$index[,2][shapefile$shx$index[,2]==2]<-6
	#Calculate new offset values for each record from beginning of shp file
	new.offset <- cumsum(c(50,shapefile$shx$index[,2]))
	#Replace the existing offsets with the new ones
	shapefile$shx$index[,1] <- new.offset[1:(length(new.offset)-1)]
	#Set the content.lenghts back to their respective values
	shapefile$shx$index[,2][shapefile$shx$index[,2]==14]<-10
	shapefile$shx$index[,2][shapefile$shx$index[,2]==6]<-2

	shapefile
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Write out SHX file (index)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write.shx <- function(shx, out.name) {
	
	outfile<-file(out.name,"wb")
	
	#Header
	writeBin(as.integer(9994), outfile, endian="big")
	writeBin(as.integer(rep(0,5)), outfile, endian="big")
	writeBin(as.integer(shx$header$file.length), outfile, endian="big")
	writeBin(as.integer(1000), outfile, endian="little")
	writeBin(as.integer(shx$header$shape.type), outfile, endian="little")
	writeBin(as.double(shx$header$xmin), outfile, endian="little")
	writeBin(as.double(shx$header$ymin), outfile, endian="little")
	writeBin(as.double(shx$header$xmax), outfile, endian="little")
	writeBin(as.double(shx$header$ymax), outfile, endian="little")
	writeBin(as.double(shx$header$zmin), outfile, endian="little")
	writeBin(as.double(shx$header$zmax), outfile, endian="little")
	writeBin(as.double(shx$header$mmin), outfile, endian="little")
	writeBin(as.double(shx$header$mmax), outfile, endian="little")
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (record in 1:nrow(shx$index)) {
		#Record Offset
		writeBin(as.integer(shx$index[record,1]), outfile, endian="big")
		#Content Length
		writeBin(as.integer(shx$index[record,2]), outfile, endian="big")
	}
	close(outfile)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Read DBF format
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

read.dbf <- function(dbf.name, header=FALSE) {
	
	#Load foreign package for read.dbf
	library(foreign)
	
	#Read header if needed
	if(header) {
		infile<-file(dbf.name,"rb")

		#Header
		file.version <- readBin(infile,integer(), 1, size=1, endian="little")
		file.year <- readBin(infile,integer(), 1, size=1, endian="little")
		file.month <- readBin(infile,integer(), 1, size=1, endian="little")
		file.day <- readBin(infile,integer(), 1, size=1, endian="little")
		num.records <- readBin(infile,integer(), 1, size=4, endian="little")
		header.length <- readBin(infile,integer(), 1, size=2, endian="little")
		record.length <- readBin(infile,integer(), 1, size=2, endian="little")
		file.temp <- readBin(infile,integer(), 20, size=1, endian="little")
		header <- list(file.version,file.year, file.month, file.day, num.records, header.length, record.length)
		names(header) <- c("file.version","file.year","file.month","file.day","num.records","header.length","record.length")
		rm(file.version,file.year, file.month, file.day, num.records, header.length, record.length)

		#Calculate the number of fields
		num.fields <- (header$header.length-32-1)/32
		field.name <- NULL
		field.type <- NULL
		field.length <- NULL
		field.decimal <- NULL

		#Field Descriptions (32 bytes each)
		for (i in 1:num.fields) {
			field.name.test <- readBin(infile,character(), 1, size=10, endian="little")
			field.name <- c(field.name,field.name.test)
			if (nchar(field.name.test)!=10) {
				file.temp <- readBin(infile,integer(), 10-(nchar(field.name.test)), 1, endian="little")
			}	
			field.type <- c(field.type,readChar(infile, 1))
			file.temp <- readBin(infile,integer(), 4, 1, endian="little")
			field.length <- c(field.length,readBin(infile,integer(), 1, 1, endian="little"))
			field.decimal <- c(field.decimal, readBin(infile,integer(), 1, 1, endian="little"))
			file.temp <- readBin(infile,integer(), 14, 1, endian="little")
		}

		#Create a table of the field info
		fields <- data.frame(NAME=field.name,TYPE=field.type,LENGTH=field.length,DECIMAL=field.decimal)
		#Set all fields with length<0 equal to correct number of characters
		fields$LENGTH[fields$LENGTH<0]<-(256+fields$LENGTH[fields$LENGTH<0])
		#Read in end of attribute descriptions terminator - should be integer value 13
		file.temp <- readBin(infile,integer(), 1, 1, endian="little")
		#Increase the length of field 1 by one to account for the space at the beginning of each record	
		fields$LENGTH[1]<-fields$LENGTH[1]+1
		#Add fields to the header list
		header <- c(header,fields=NULL)
		header$fields <- fields

		#Close file connection
		close(infile)
	}
	
	#Read in dbf
	dbf <- get("read.dbf","package:foreign")(dbf.name)

	#Return the dbf as a list with a data.frame and a header list
	list(dbf=dbf, header=header)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Write out DBF format
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.dbf <- function(dbf, out.name, arcgis=FALSE) {
	
	#Load foreign package for write.dbf
	library(foreign)
	
	#If output intended for ArcGIS, replace "." with "_" for column names in header
	if (arcgis==T) {
		colnames(dbf$dbf) <- gsub("\\.","_",colnames(dbf$dbf))
	}
	
	#Write out dbf
	get("write.dbf","package:foreign")(dbf$dbf, out.name)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Function to Read in a Shapefile
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

read.shapefile <- function(shape.name) {
	shp.data <- read.shp(paste(shape.name, ".shp", sep=""))
	shx.data <- read.shx(paste(shape.name, ".shx", sep=""))
	dbf.data <- get("read.dbf","package:shapefiles")(paste(shape.name, ".dbf", sep=""))
	shapefile <- list(shp=shp.data,shx=shx.data,dbf=dbf.data)		
	shapefile
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Function to Write out a Shapefile
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.shapefile <- function(shapefile, out.name, arcgis=FALSE) {
	write.shp(shapefile$shp, paste(out.name, ".shp", sep=""))
	write.shx(shapefile$shx, paste(out.name, ".shx", sep=""))
	get("write.dbf","package:shapefiles")(shapefile$dbf, paste(out.name, ".dbf", sep=""), arcgis)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Function to Add the X and Y Coordinates to the DBF of a shapefile
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

add.xy <- function(shapefile) {
	if (shapefile$shp$header$shape.type !=1) {
		stop ("Must be point shapefile")
	}
	shapefile$dbf[[1]]$XCOORD <- shapefile$shp[[1]][,2]
	shapefile$dbf[[1]]$YCOORD <- shapefile$shp[[1]][,3]
	shapefile
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Function to scale the X and Y Coordinates of the shapefile
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

scaleXY <- function(shapefile, scale.factor) {
	#Scale the bounding box of the shapefile by the scale factor
	shapefile$shp$header$xmin <- shapefile$shp$header$xmin/scale.factor
	shapefile$shp$header$ymin <- shapefile$shp$header$ymin/scale.factor
	shapefile$shp$header$xmax <- shapefile$shp$header$xmax/scale.factor
	shapefile$shp$header$ymax <- shapefile$shp$header$ymax/scale.factor
	
	#Scale the bounding box of the shx file by the scale factor
	shapefile$shx$header$xmin <- shapefile$shx$header$xmin/scale.factor
	shapefile$shx$header$ymin <- shapefile$shx$header$ymin/scale.factor
	shapefile$shx$header$xmax <- shapefile$shx$header$xmax/scale.factor
	shapefile$shx$header$ymax <- shapefile$shx$header$ymax/scale.factor
	
	#Scale the X and Y if a point shape
	if (shapefile$shp$header$shape.type == 1) {
		shapefile$shp[[1]][,2] <- shapefile$shp[[1]][,2]/scale.factor
		shapefile$shp[[1]][,3] <- shapefile$shp[[1]][,3]/scale.factor
	}
	#Scale the X and Y point values and bounding box if a line or polyon shape
	if (shapefile$shp$header$shape.type==3 || shapefile$shp$header$shape.type==5) {
		for (shape in 1:length(shapefile$shp$shp)) {
			shapefile$shp[[1]][[shape]]$points[,1] <- shapefile$shp[[1]][[shape]]$points[,1]/scale.factor
			shapefile$shp[[1]][[shape]]$points[,2] <- shapefile$shp[[1]][[shape]]$points[,2]/scale.factor
			shapefile$shp[[1]][[shape]]$box <- shapefile$shp[[1]][[shape]]$box/scale.factor
		}
	}
	shapefile
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Converts a simple point, polyLine, or polygon data frame into a shapefile. 
#The shpTable data frame must have three columns in this order: Id, X, and Y. 
#The attTable data frame must have a field to join on - identified by field
#The type argument is either 1 (point), 3 (polyline) or 5 (polygon).
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
convert.to.shapefile <- function(shpTable, attTable, field, type) {
	
	#Check to ensure field in attTable
		if(!(field %in% colnames(attTable))) {
			stop("Field not in attTable column names")
	}
	
	#Check to ensure the same Id set in shpTable and attTable
	if(length(unique(shpTable[,1])) != length(attTable[,field])) {
		stop("Different number of unique Ids in shpTable versus attTable")
	}
	
	if(any(sort(unique(shpTable[,1])) != sort(attTable[,field]))) {
		stop("Id set in shpTable not the same as Id set in attTable") 
	}
	
	#Build dbf - no longer uses header since "foreign.write.dbf" does not
	dbf <- list(dbf=attTable,header=NULL)

	#Build shp file

	#shp list
	shpTable[,c(2,3)] <- round(shpTable[,c(2,3)], 5)
	colnames(shpTable) <- c("Id","X","Y")
	v.content.length <- vector()

	#point
	if(type == 1) {
	
		#loop for each entry in the dbf and build a shape
		for (i in 1:nrow(attTable)) {
		
			#Get shape entries for attTable Id i
			shpTable.i <- shpTable[which(shpTable[,1] == attTable[i,field]), -1] 
			p.content.length <- 10 #shape.type and 2 doubles
			v.content.length[i] <- p.content.length
			if(i == 1) {
				shp <- data.frame(
					record = i,
					x = shpTable.i[,1],
					y = shpTable.i[,2],
					shape.type = type)
			} else {
				shp2 <- data.frame(
					record = i,
					x = shpTable.i[,1],
					y = shpTable.i[,2],
					shape.type = type)
				shp <- rbind(shp, shp2)
				rm(shp2)
			}
		}
	}	
	
	#polyLine and polygon shp
	if(type == 3 || type == 5) {
		
		#shp file
		shp<-list()
		
		#loop for each entry in the dbf and build a shape
		for (i in 1:nrow(attTable)) {
		
			#Get shape entries for attTable Id i
			shpTable.i <- shpTable[which(shpTable[,1] == attTable[i,field]), -1] 
			num.points <- nrow(shpTable.i)
			v.box <- c(min(shpTable.i[,1]),min(shpTable.i[,2]),max(shpTable.i[,1]),max(shpTable.i[,2]))
			names(v.box) <- c("xmin","ymin","xmax","ymax")
			p.content.length <- 24 + (num.points * 8)
			v.content.length[i] <- p.content.length
			shp[[i]]<-list(record = i,
				content.length = p.content.length,
				shape.type = type,
				box = v.box,
				num.parts = 1, #Only 1 part
				num.points = nrow(shpTable.i),
				parts = 0, 
				points = shpTable.i)
		}
	}

	#Collect up all shp info
	header<-list(file.code = 9994,
		file.length = 50 + sum(v.content.length + 4),
		file.version = 1000,
		shape.type = type,
		xmin = min(shpTable[,2]),
		ymin = min(shpTable[,3]),
		xmax = max(shpTable[,2]),
		ymax = max(shpTable[,3]),
		zmin = 0,
		zmax = 0,
		mmin = 0,
		mmax = 0)
	shp <- list(shp=shp, header=header)

	#Build shx file

	#shx list
	shx<-list(index = cbind(Offset=c(50 ,(cumsum(v.content.length + 4) + 50)[-length(v.content.length)]),
		Length = v.content.length),
		header = shp$header)

	shx$header$file.length <- 50 + 4 * nrow(shx$index)

#Return shapefile
shp.file <- list(shp=shp, shx=shx, dbf=dbf)
shp.file

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Convert read.shp result to simple data.frame format (Id, X, Y)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
convert.to.simple <- function(shp) {
	
	#Determine type of shp
	if(shp$header$shape.type == 1) {
		#Get points
		allPoints <- shp$shp[,c("record","x","y")]
		colnames(allPoints) <- c("Id", "X", "Y")
	
	} else {
		#Get points
		allPoints <- do.call(rbind, lapply(shp$shp, function(x) x$points[,c("X","Y")]))
		
		#Calculate index and add to points
		indexes <- 1:length(shp$shp)
		repNum <- lapply(shp$shp, function(x) nrow(x$points))
		allPoints <- cbind(Id=rep(indexes, repNum), allPoints)
	}
	
	#Return allPoints
	allPoints
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Change Index to a Field from the DBF
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Changes the Index field of the shpTable to the new vector of fields
#Will replace the records in order, so Index 1 = element 1 in the new vector
change.id <- function(shpTable, newFieldAsVector) {
	
	#Check to ensure the same Id set in shpTable and attTable
	if(length(unique(shpTable[,1])) != length(newFieldAsVector)) {
		stop("Different number of unique Ids in shpTable versus new field")
	}

	#Create new vector
	count <- table(shpTable[,1])
	shpTable[,1] <- rep(newFieldAsVector, count)
	
	#Return result
	shpTable
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Douglas-Peucker Polyline Simplification
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Ben Stabler, bstabler@ptvamerica.com, March 2005
#Douglas, D. and Peucker, T. (1973). "Algorithms for 
#the reduction of the number of points required to 
#represent a digitized line or its caricature." 
#The Canadian Cartographer 10(2). 112-122.

#Currently uses the line, not the line segment to 
#determine the distance of the points from the line.
#See http://www.lgc.com/resources/Doug_Peucker.pdf
#for more information.  This can result in the 
#omission of extreme "outlier-like" points. Try:
#points<-list(x=c(200,0,400),y=c(0,20,0))
#plot(points, type="l")
#lines(dp(points, 50), type="l", col="blue")

#Simple Example
#x <- c(5,3,4,1,8,9,10,11)
#y <- c(6,4,2,1,1,5,2,3)
#points <- list(x=x,y=y)
#plot(points, type="l")
#lines(dp(points, 2), type="l", col="blue")

####################################################

dp <- function(points, tolerance) {

 #Convert to lowercase
 names(points) <- tolower(names(points))

 #Calculate distance between two points
 distance <- function(x1,x2,y1,y2) {
 	sqrt((x2-x1)^2 + (y2-y1)^2)
 }
 
 #Calculate equation of a line from two points
 equationOfLine <- function(x1,x2,y1,y2) {
 	slope <- (y2-y1)/(x2-x1)
 	b <- y1 - slope * x1
 	c(slope,b)
 }
 
 #Calculate y-intercept from a point and a slope
 calcB <- function(slope, px, py) {
 	py - (slope * px) 
 }
 
 #Calculate intercept of two lines from their slope and y-intercept
 intercept <- function(s1,b1,s2,b2) {
 	x1 <- (b2-b1)/(s1-s2)
 	y1 <- s1 * x1 + b1
 	c(x1,y1)
 }

 #Setup vector to mark points to keep
 keep <- rep(F, length(points$x))
 keep[1] <- T
 keep[length(keep)] <- T

 #Function definition to simplify points
 simplify <- function(start, end, tol=tolerance) {

  #Calculate intermediate point distances 
  if (length(points$x[start:end]) > 2) {

 	#Avoid Inf slope
 	if( points$x[start] ==  points$x[end] ) { points$x[start] <- points$x[start] - 0.0000001 }
 	if( points$y[start] ==  points$y[end]) { points$y[start] <- points$y[start] - 0.0000001 }
 
	#Calculate equation of line of middle points
 	line <- equationOfLine( points$x[start], points$x[end], points$y[start], points$y[end])
	
	#Calculate y-intercepts
 	b <- mapply(function(x,y) calcB(-1/line[1], x, y), points$x[(start+1):(end-1)], points$y[(start+1):(end-1)])
	
	#Calculate intercepts with with start-end line
 	ints <- sapply(b, function(x) intercept(line[1], line[2], -1/line[1], x), simplify=F)

	#Calculate distances of points from line
 	distances <- mapply(function(x,y,z) distance(x[[1]], y, x[[2]], z), ints, 
		points$x[(start+1):(end-1)], points$y[(start+1):(end-1)])

	#If any point greater than tolerance split at max distance point and apply to each side
 	if (any(distances >= tol)) {
		keep[which.max(distances)+start] <<- T
		#print(which.max(distances)+start)
		
		simplify(start, which.max(distances)+start)
		simplify(which.max(distances)+start, end)
	}

  }

 }

 #Start simplification 
 simplify(1, length(points$x))

 #Return simplified points
 list(x=points$x[keep], y=points$y[keep])
}

