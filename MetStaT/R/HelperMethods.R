# 		HelperMethods: selection of methods to help tools in the 'MetStaT' package  
#		Copyright (C) 2012 Tim Dorscheidt (see ScalePip, mldivide and mrdivide methods for additional copyright details)
#		
#		This program is free software: you can redistribute it and/or modify
#		it under the terms of the GNU General Public License as published by
#		the Free Software Foundation, either version 3 of the License, or
#		(at your option) any later version.
#		
#		This program is distributed in the hope that it will be useful,
#		but WITHOUT ANY WARRANTY; without even the implied warranty of
#		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#		GNU General Public License for more details.
#		
#		You should have received a copy of the GNU General Public License
#		along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# 		Email: g.zwanenburg@uva.nl ('MetStaT' contact person) or tdorscheidt@gmail.com
###############################################################################

# This method is an adjusted version of R's 'scale' method in the 'base' package, GPL license
MetStaT.ScalePip <- function(x.input, center = TRUE, scale = TRUE, quietly = FALSE)
{
	options(warn=-1)
	no.col.x.input <- ncol(x.input)
	if (is.null(no.col.x.input)) {no.col.x.input <- 1}
	tryCatch({
				x <- matrix(as.numeric(x.input),ncol=no.col.x.input)
			}, error = function(ex) {
				bad.matrix <- x.input;stop(ex);
			})
	colnames(x) <- colnames(x.input)
	rownames(x) <- rownames(x.input)
	options(warn=0)
	
	x.scaled <- list()
	nc <- ncol(x)
	if (is.null(center)) center <- FALSE
	if (is.character(center) && center=="true") center <- TRUE
	if (is.character(center) && center=="false") center <- FALSE
	if (is.character(scale) && scale=="true") scale <- TRUE
	if (is.character(scale) && scale=="false") scale <- FALSE
	
	center.description <- center
	if (is.logical(center)) {
		if (center) {
			center.description <- "Around mean. "
			center <- colMeans(x, na.rm=TRUE)
			x <- sweep(x, 2L, center, check.margin=FALSE)
		}
		else {
			x.scaled$description <- paste(x.scaled$description,"Not centered. ",sep="")
			not.centered <- matrix(rep(0, nc),nrow=1)
			colnames(not.centered) <- colnames(x)
			x.scaled$center.vector <- not.centered
		}
	}	else if (is.numeric(center) && (length(center) == nc)) {
		center.description <- "Manual input by user used. "
		x <- sweep(x, 2L, center, check.margin=FALSE)
	} else {
		stop("length of 'center' must equal the number of columns of 'x'")
	}
	
	if(is.numeric(center)) {
		x.scaled$description <- paste(x.scaled$description,"Centered: ",center.description,sep="")
		center <- matrix(center,nrow=1)
		colnames(center) <- colnames(x)
		x.scaled$center.vector <- center
	}
	
	if (is.null(scale)) scale <- FALSE
	if (is.logical(scale)) {
		if (scale) {
			scale="stdev"
		}
	}
	
	scale.description <- scale
	if (is.logical(scale)) {
		x.scaled$description <- paste(x.scaled$description,"Not scaled. ",sep="")
		not.scaled <- matrix(rep(1, nc),nrow=1)
		colnames(not.scaled) <- colnames(x)
		x.scaled$scale.vector <- not.scaled
	} else if (is.character(scale)) {
		scale <- tolower(scale)
		if (scale == "stdev" || scale=="auto") {
			f <- function(v) {
				v <- v[!is.na(v)]
				sqrt(sum(v^2) / max(1, length(v) - 1L))
			}
		}	else if (scale == "pareto") {
			f <- function(v) {
				v <- v[!is.na(v)]
				sqrt(sqrt(sum(v^2) / max(1, length(v) - 1L)))
			}
		}
		scale <- apply(x, 2L, f)
		x <- sweep(x, 2L, scale, "/", check.margin=FALSE)
	}	else if (is.numeric(scale) && length(scale) == nc) {
		scale.description <- "Manual input by user used."
		x <- sweep(x, 2L, scale, "/", check.margin=FALSE)
	}	else {
		stop("length of 'scale' must equal the number of columns of 'x'")
	}
	
	if(is.numeric(scale)) {
		x.scaled$description <- paste(x.scaled$description,"Scaled: ",scale.description,".",sep="")
		scale <- matrix(scale,nrow=1)
		colnames(scale) <- colnames(x)
		x.scaled$scale.vector <- scale
	}
	
	x.scaled$data <- x
	if (!quietly) {print(x.scaled$description)}
	x.scaled
}

# method that ensures all plot devices are killed
MetStaT.KillAllDevices <- function() {
	no.devices <- length(dev.list())
	for (d in 1:no.devices) {
		dev.off(which=dev.list()[1])
	}
}

# Method that converts values to numerical classes
MetStaT.ConvertToNumericClasses <- function(data, cols = NULL, new.classes = NULL) {
	ModNoZero <- function(a,b) { # for this function I need an adjusted modulo conversion
		res <- a%%b
		if (res==0) {res <- b}
		res
	}
	data.dim <- dim(data)
	if (is.null(data.dim)) {data <- as.matrix(data); data.dim <- dim(data)}
	if (is.null(cols)) {
		cols = c(1:dim(data)[2])
	} else if (!is.numeric(cols)) {
		eval(parse(text = paste(sep="","cols <- c(",cols,")")))
	}
	data.out <- matrix(NA,nrow=dim(data)[1],ncol=dim(data)[2])
	if (!is.null(new.classes)) {
		if (is.numeric(new.classes)) {
			class.labels <- new.classes
		} else {
			eval(parse(text = paste(sep="","class.labels <- c(",new.classes,")")))
		}
	}
	for (c in cols) {
		if (is.null(new.classes) || length(new.classes)==0) { # user did not not define new classes, create default ones
			no.classes.needed <- length(unique(data[,c])) # count unique values in current data-column
			if (no.classes.needed==2) {
				class.labels <- c(-1,1);
			} else {class.labels <- c(1:no.classes.needed);}
		}
		unique.values.present <- unique(data[,c])
		class.index <- 1
		for (u in unique.values.present) {
			data.out[data[,c]==u,c] <- class.labels[ModNoZero(class.index,length(class.labels))]
			class.index <- class.index + 1
		}
	}
	data.out
}

# Method that reads a data file and outputs a matrix, whereby the first line of the file is interpreted as column names
MetStaT.ReadFileToHeaderMatrix <- function(file.to.read, file.contains.header = TRUE, file.contains.row.names = FALSE, rows = "", cols = "", separator = "\t", force.numeric = FALSE) {
	result <- as.matrix(read.delim(file = file.to.read, header = file.contains.header, fill = FALSE, sep = separator),byrow = T)
	if (file.contains.row.names) {
		row.names <- result[,1,drop=FALSE]
	}
	if (rows!="" || cols!="") {
		if (rows!="") rows = paste(sep="","c(",rows,")")
		if (cols!="") cols = paste(sep="","c(",cols,")")
		cutexpression <- paste(sep="","result <- result[",rows,",",cols,",drop=FALSE]")
		eval(parse(text=cutexpression))
	}
	if (file.contains.header == FALSE) {
		nc <- ncol(result)
		colname = paste(sep="","c",1:nc)
		colnames(result) <- colname
	}
	if (file.contains.row.names) {
		rownames(result) <- row.names
	}
	if (mode(result)!="numeric" && force.numeric) {
		options(warn=-1)
		if (try( sum(is.na(as.numeric(result))) ) != 0) {
			result.numeric <- matrix(nrow = dim(result)[1], ncol = dim(result)[2])
			for (col in 1:dim(result)[2]) {
				if (try( sum(is.na(as.numeric(result[,col]))) ) != 0) {
					result.numeric[,col] <- NA
				} else {
					result.numeric[,col] <- as.numeric(result[,col])
					result.numeric[,col] <- as.numeric(result.numeric[,col])
				}
			}
			colnames(result.numeric) <- colnames(result)
			rownames(result.numeric) <- rownames(result)
			result <- result.numeric      
		}
		else {
			result.numeric <- matrix(as.numeric(result),nrow = nrow(result), ncol = ncol(result))
			colnames(result.numeric) <- colnames(result)
			rownames(result.numeric) <- rownames(result)
			result <- result.numeric
		}
		options(warn=0)
	}
	result
}

# Convert a matrix to numerical values (or NA) without losing any meta-data
MetStaT.ConvertToNumeric <- function(matrix) {
	if (is.numeric(matrix)) {return(matrix)}
	output <- matrix(nrow=dim(matrix)[1], ncol=dim(matrix)[2])
	names(output) <- names(matrix)
	forced.na <- FALSE
	options(warn=-1)
	for (row in 1:dim(matrix)[1]) {
		for (col in 1:dim(matrix)[2]) {
			if (is.na(matrix[row,col]) || toupper(matrix[row,col])=="NA") {
				output[row,col]=NA
			} else if (is.na(as.numeric(matrix[row,col]))) {
				output[row,col]=NA; forced.na <- TRUE
			} else {
				output[row,col]=as.numeric(matrix[row,col])
			}
		}
	}
	options(warn=0)
	if (forced.na) {
		print("Caution: your data contains non-numeric values. Certain tools may act unpredictable. Perhaps you ticked the incorrect option for: 'First row in dataset contains column names:'.")
	}
	output
}

# Write a matrix with column names to a value separated file
MetStaT.CreateFileFromHeaderMatrix <- function(file.to.create, header.matrix, rownames = FALSE) {
	write.table(header.matrix,file.to.create, sep = "\t", col.names = TRUE, row.names = rownames)
}

# Write a matrix with column and row names to a value separated file
MetStaT.CreateFileFromHeaderRowMatrix <- function(file.to.create, row.header.matrix, description="") {
	row.header.matrix<-as.table(row.header.matrix)
	if (is.null(rownames(row.header.matrix))) {
		rownames(row.header.matrix) <- paste(sep="","r",1:dim(row.header.matrix)[1])
	}
	if (is.null(colnames(row.header.matrix))) {
		colnames(row.header.matrix) <- paste(sep="","c",1:dim(row.header.matrix)[2])
	}
	clean.matrix <- rbind(colnames(row.header.matrix),row.header.matrix)
	clean.matrix <- cbind(c(description,rownames(row.header.matrix)),clean.matrix)
	write.table(clean.matrix,file.to.create, sep = "\t", col.names = FALSE, row.names = FALSE)
	clean.matrix
}

# Write all graphical plots to a file that are the result of certain R expressions
MetStaT.PlotToFile <- function(zipfile.name, filename.no.ext, plot.expressions, file.type, ...) {
	if (file.type=="pdf") {
		filename.with.ext <- paste(sep="",filename.no.ext,".pdf")
		pdf(filename.with.ext, ...) # width = output.width, height = output.height, ...)
	} else if (file.type=="bmp") {
		filename.with.ext <- paste(sep="",filename.no.ext,".bmp")
		bmp(filename.with.ext, ...) #  width = output.width, height = output.height, ...)
	} else if (file.type=="jpeg") {
		filename.with.ext <- paste(sep="",filename.no.ext,".jpeg")
		jpeg(filename.with.ext, ...) #  width = output.width, height = output.height, ...)
	} else if (file.type=="png") {
		filename.with.ext <- paste(sep="",filename.no.ext,".png")
		png(filename.with.ext, ...) #  width = output.width, height = output.height, ...)
	} else if (file.type=="tiff") {
		filename.with.ext <- paste(sep="",filename.no.ext,".tiff")
		tiff(filename.with.ext, ...) #  width = outputW.width, height = output.height, ...)
	}
	eval(parse(text=plot.expressions))
	dev.off(which=dev.cur())
	filename.with.ext
	system(paste(sep=" ","zip", zipfile.name, filename.with.ext))
	system(paste(sep=" ","rm", filename.with.ext))
}

# Write some common R objects to a human readable file
MetStaT.WriteDataObjectToFile <- function(filename, data) {
	GetCharsToWriteInternalFunction <- function(data) {  
		if (is.list(data)) {
			chars.to.write <- c()
			names.of.children <- names(data)
			if (is.null(names.of.children)) {names.of.children <- 1:length(data)}
			for (child.index in 1:length(data)) {
				chars.to.write <- c(chars.to.write, names.of.children[child.index])
				chars.to.write <- c(chars.to.write, GetCharsToWriteInternalFunction(data[[child.index]]))
			}
		} else if (is.matrix(data)) {
			chars.to.write <- c()
			if (!is.null(colnames(data))) {
				if (!is.null(rownames(data))) {
					chars.to.write <- c(chars.to.write, paste(" ", paste(colnames(data),collapse="\t"), sep="\t"))
				} else {
					chars.to.write <- c(chars.to.write, paste(colnames(data),collapse="\t"))
				}
			}
			if (!is.null(rownames(data))) {
				chars.to.write <- c(chars.to.write, paste(rownames(data),apply(data,1,paste,collapse="\t"),sep="\t"))
			} else { chars.to.write <- c(chars.to.write, apply(data,1,paste,collapse="\t")) }
		} else if (is.numeric(data)) {
			chars.to.write <- c()
			if (!is.null(colnames(data))) {chars.to.write <- c(chars.to.write, paste(colnames(data),collapse="\t"))}
			chars.to.write <- c(chars.to.write,paste(data,collapse="\t"))
		} else if (is.character(data)) {
			chars.to.write <- data
		} else {chars.to.write <- "Data type unknown"}
		chars.to.write
	}
	
	file.conn<-file(filename)
	writeLines(GetCharsToWriteInternalFunction(data), file.conn)
	close(file.conn)
}

# convert a string with numerical values and comma's or ranges of numerical values to a vector
MetStaT.ConcatWithStringPars <- function(string) {
	values = eval(parse(text=paste(sep="","c(",string,")")))
	values
}

# export the results of certain R expressions to data-files, and wrap all results in a zip file
MetStaT.ExportDataToFile <- function(zipfile.name, filename.no.ext, data.expressions, file.type) {
	for (i in 1:length(data.expressions)) {
		if (length(data.expressions)==1) {
			filename.with.ext <- paste(sep="",filename.no.ext,".",file.type)
		} else {filename.with.ext <- paste(sep="",filename.no.ext,"(",i,").",file.type)}
		data.result <- eval(parse(text=data.expressions[i]))
		MetStaT.WriteDataObjectToFile(filename.with.ext,data.result)
		system(paste(sep=" ","zip", zipfile.name, filename.with.ext))
		system(paste(sep=" ","rm", filename.with.ext))
	}
}

# function to remove columns that contain one or more NAs
MetStaT.RemoveNaColumns <- function(input.matrix, rows.to.ignore = NULL) {
	na.columns <- c()
	matrix.to.scan <- input.matrix
	if (!is.null(rows.to.ignore)) {matrix.to.scan <- matrix.to.scan[-rows.to.ignore,,drop=FALSE]}
	for (c in 1:dim(matrix.to.scan)[2]) { 
		if (any(is.na(matrix.to.scan[,c]))) {na.columns <- c(na.columns,c)}
	}
	input.matrix[,-na.columns,drop=FALSE]
}

# the following function trims string that are too long, and ends them with a user-defined trim ending and possibly some of the ending of the input-string
MetStaT.TrimCustom <- function(text, max.length=5, trim.ending="..", include.ending.length=0) {
	if ((nchar(trim.ending)+include.ending.length)>=max.length) {stop("Trimmed output has ending which is longer than total allowed length of output: length of trimEnding and includeEndingLength exceeds maxLength.")} 
	if (nchar(text)>max.length) {
		if (include.ending.length>0) {
			ending.text <- substr(text,nchar(text)-include.ending.length+1,nchar(text))
		}
		text <- paste(substr(text,1,max.length-(nchar(trim.ending)+include.ending.length)),trim.ending,ending.text,sep="")
		
	}  
	text
}

# helper-function that calculates how often a list of pre-defined class-types occurs in an array
MetStaT.GetFreqTable <- function(classes.to.check, class.types = NULL) {
	if (is.null(class.types)) {class.types <- unique(classes.to.check)}
	result <- array(NA,length(class.types))
	names(result) <- class.types
	for (i in 1:length(class.types)) {
		result[i] <- sum(classes.to.check==class.types[i])
	}
	result
}

# function from 'pracma' package by Hans W Borchers, GPL license
MetStaT.mldivide <- function(A, B) { # backslash operator \
	stopifnot(is.numeric(A) || is.complex(A),
			is.numeric(B) || is.complex(B))
	if (is.vector(A)) A <- as.matrix(A)
	if (is.vector(B)) B <- as.matrix(B)
	if (nrow(A) != nrow(A))
		stop("Matrices 'A' and 'B' must have the same number of rows.")
	
	qr.solve(A, B)
}

# function from 'pracma' package by Hans W Borchers, GPL license
MetStaT.mrdivide <- function(A, B) { # forward slash operator /
	stopifnot(is.numeric(A) || is.complex(A),
			is.numeric(B) || is.complex(B))
	if (is.vector(A)) A <- t(A)
	if (is.vector(B)) B <- t(B)
	if (ncol(A) != ncol(A))
		stop("Matrices 'A' and 'B' must have the same number of columns.")
	
	t(MetStaT.mldivide(t(B), t(A)))
}

# Method used to find all possible pairings of principal components within a set 
MetStaT.GetPcTuples <- function(no.pc) {
	pcs <- 1:no.pc
	list.of.tuples <- list()
	index.pc1 <- 1
	while(index.pc1 <= no.pc) {
		index.pc2 <- index.pc1+1  
		while(index.pc2 <= no.pc) {
			list.of.tuples <- c(list.of.tuples,list(c(index.pc1,index.pc2)))    
			index.pc2 <- index.pc2 + 1 
		}
		index.pc1 <- index.pc1 + 1  
	}
	list.of.tuples
}
