library(methods);

### split set of numbers into increasing sequences

.index.splitter = function(i) {
	parts = which(c(TRUE, diff(i)!=1,TRUE));
	ind = list( start = i[parts[1:(length(parts)-1)]], len = diff(parts), ix = parts, n = length(parts)-1);
	return(ind);
}
# Example
# .index.splitter(c(1,2,3, 3,4,5, -3,-2,-1))
# Returns:
# 
# $start - first value in each sequence,      1 3 -3
# $len - length of each sequence,             3 3 3
# $ix - start position in the input sequence, 1 4 7 10 
# $n - number of sequences,                   3

### Reference class for the filematrix
setRefClass("filematrix",
	fields = list( 
		fid = "list",                # list with file handle for the binary data file
		nr = "numeric",              # number of rows, access via nrow(fm) and dim(fm)
		nc = "numeric",              # number of columns, access via ncol(fm) and dim(fm)
		type = "character",          # type of values in the matrix (double,integer,logical,raw)
		size = "integer",            # size of each matrix value in the file (1,2,4,8)
		caster = "function",         # function transforming data into the matrix data type
		info.filename = "character", # file name for file matrix description
		rnames = "character",        # row names, access via rownames(fm)
		cnames = "character",        # column names, access via colnames(fm)
		rnamefile = "character",     # file with row names
		cnamefile = "character"      # file with column names
	),
	methods = list(
		# Set the caster function based on data "type"
		setCaster = function() {
			caster <<- switch(type,
				double = as.double,
				integer = as.integer,
				logical = as.logical,
				raw = as.raw,
				stop("Unknown data type: ",type));
		},
		# Initialize all variables in the class. Called automatically upon creation.
		initialize = function() {
			.self$fid = list();
			.self$nr = 0L;
			.self$nc = 1L;
			.self$type = "double";
			.self$size = 8L;
			.self$setCaster();
			.self$info.filename = "";
			.self$rnames = character();
			.self$cnames = character();
			.self$rnamefile = "";
			.self$cnamefile = "";
		},
		# Close the file matrix object. Access via close(fm)
		close = function() {
			if( length(.self$fid)>0 ) {
				base::close.connection( .self$fid[[1]] );
				.self$fid = list();
			} else {
				warning("Inactive filematrix object.");
			}
			initialize();
		},
		# Is open?
		isOpen = function() {
			return( length(.self$fid)>0 );
		},
		# This method is called when the object is being printed.
		show = function() {
			if(length(fid)>0) {
				cat(nr, "x", nc, "filematrix object", "\n");
			} else {
				cat("Inactive filematrix object","\n");
			}
		},
		# methods for reading from and writing to the descriptor file "info.filename"
		loadInfo = function() {
			info = readLines( .self$info.filename);
			keep = grep(x=info, pattern="=", fixed=TRUE);
			info = info[keep];
			ss = strsplit(info, split="=");
			lst = lapply(ss, "[[", 2);
			names(lst) = sapply(ss, "[[", 1);
			
			.self$nr = round(as.numeric(lst$nrow));
			.self$nc = round(as.numeric(lst$ncol));
			.self$size = as.integer(lst$size);
			.self$type = lst$type;
			if(!(type %in% c("double","integer","logical","raw"))) {
				stop("\"type\" must be either \"double\",\"integer\",\"logical\", or \"raw\"");
			}
			.self$setCaster();
		},
		saveInfo = function() {
			writeLines(con=.self$info.filename, text=paste0(
				"# Information file for R filematrix object","\n",
				"nrow=", .self$nr,   "\n",
				"ncol=", .self$nc,   "\n",
				"type=", .self$type, "\n",
				"size=", .self$size, "\n"));
		},
		# Get and set dimension names. Access via rownames(fm), colnames(fm), and dimnames(fm).
		getrownames = function() {
			if( length(rnames) > 0 ) {
				rn = rnames;
			} else {
				if(file.exists(rnamefile)) {
					rn = readLines(rnamefile);
					.self$rnames = rn;
				} else {
					rn = NULL;
				}
			}
			return(rn);
		},
		getcolnames = function() {
			if( length(cnames) > 0 ) {
				cn = cnames;
			} else {
				if(file.exists(cnamefile)) {
					cn = readLines(cnamefile);
					.self$cnames = cn;
				} else {
					cn = NULL;
				}
			}
			return(cn);
		},
		getdimnames = function() {
			return(list(getrownames(),getcolnames()));
		},
		setrownames = function(rn) {
			if(length(rn)>0) {
				.self$rnames = rn;
				writeLines( text=rnames, con=rnamefile);
			} else {
				if(file.exists(rnamefile))
					file.remove(rnamefile);
				.self$rnames = character();
			}
			return(invisible(.self));
		},
		setcolnames = function(cn) {
			if(length(cn)>0) {
				.self$cnames = cn;
				writeLines( text=cnames, con=cnamefile);
			} else {
				if(file.exists(cnamefile))
					file.remove(cnamefile);
				.self$cnames = character();
			}
			return(invisible(.self));
		},
		setdimnames = function(nms) {
			setrownames(nms[[1]]);
			setcolnames(nms[[2]]);
			return(invisible(.self));
		},
		# File creation functions
		create = function(filenamebase, nrow = 0, ncol = 1, type="double", size=NULL) {
			
			filenamebase = gsub("\\.desc\\.txt$", "", filenamebase);
			filenamebase = gsub("\\.bmat$", "", filenamebase);
			filenamebase = normalizePath(filenamebase, mustWork=FALSE);
			.self$rnamefile =     paste0(filenamebase, ".nmsrow.txt");
			.self$cnamefile =     paste0(filenamebase, ".nmscol.txt");
			.self$info.filename = paste0(filenamebase, ".desc.txt");
			
			if( !(type %in% c("double","integer","logical","raw")) ) {
				stop("\"type\" must be either \"double\",\"integer\",\"logical\", or \"raw\"");
			}
			if(is.null(size)) {
				.self$size = switch(type,
					double = 8L,
					integer = 4L,
					logical = 1L,
					raw = 1L,
					stop("Unknown data type: ",type));
			} else {
				.self$size = as.integer(size);
			}
			
			.self$nc = round(as.double(ncol));
			.self$nr = round(as.double(nrow));
			.self$type = type;
			setCaster();
			
			saveInfo();
			
			data.file.name = paste0(filenamebase,".bmat");
			fd = file(description=data.file.name, open="w+b");
			.self$fid = list(fd);
			if(nr*nc>0)
				writeSeq(nr*nc, 0);
		},
		open = function(filenamebase, readonly = FALSE) {
			
			filenamebase = gsub("\\.desc\\.txt$", "", filenamebase);
			filenamebase = gsub("\\.bmat$", "", filenamebase);
			filenamebase = normalizePath(filenamebase, mustWork=FALSE);
			.self$rnamefile =     paste0(filenamebase, ".nmsrow.txt");
			.self$cnamefile =     paste0(filenamebase, ".nmscol.txt");
			.self$info.filename = paste0(filenamebase, ".desc.txt");
			loadInfo();

			data.file.name = paste0(filenamebase, ".bmat");
# 			stopifnot( file.info(data.file.name)$size == nr*nc*size );
			stopifnot( file.info(data.file.name)$size >= nr*nc*size );
			fd = file(description=data.file.name, open=if(readonly){"rb"}else{"r+b"});

			.self$fid = list(fd);
		},
		createFromMatrix = function(filenamebase, mat, size=NULL) {
			mat = as.matrix(mat);
			create(filenamebase=filenamebase, nrow=nrow(mat), ncol=ncol(mat), type=typeof(mat), size=size);
			setdimnames(dimnames(mat));
			writeAll(mat);
			return(invisible(.self));
		},
		# Data access routines. Access via fm[].
		# Arguments are assumed to be round (integers, possibly > 2^32).
		# No check if the object is closed.
		# Both are checked in fm[] interface.

		# fm[start:(start+len-1)]
		readSeq = function(start, len) {
			stopifnot( start>=1 );
			stopifnot( start+len-1 <= nr*nc );
			seek(con=fid[[1]], where=(start-1)*size, rw="read");
			rez = readBin(con=fid[[1]], n=len, what=type, size=size, endian="little");
			return(rez);
		},
		# fm[start:(start+length(value)-1)] = value
		writeSeq = function(start, value) {
			stopifnot( start >= 1L );
			stopifnot( start+length(value)-1 <= nr*nc );
			seek(con=fid[[1]], where=(start-1L)*size, rw="write");
			writeBin(con=fid[[1]], object=caster(value), size=size, endian="little");
			flush(fid[[1]]);
			return(invisible(.self));
		},
		# fm[,start:(start+num-1)]
		readCols = function(start, num) {
			rez = readSeq( (start-1)*nr+1, num*nr );
			if((num > 1) | (nr < 2^31))
				dim(rez) = c(nr, num);
			return(rez);
		},
		# fm[,start:(start+ncol(value)-1)] = value
		writeCols = function(start, value) {
			stopifnot( (length(value)%%nr)==0 );
			writeSeq(	(start-1)*nr+1, value);
			return(invisible(.self));
		},
		# fm[i, j:(j+num-1)]
		readSubCol = function(i, j, num) {
			stopifnot( i>=1  );
			stopifnot( j>=1  );
			stopifnot( i<=nr );
			stopifnot( j<=nc );
			stopifnot( i + num - 1 <= nr );
			rez = readSeq( (j-1)*nr+i, num);
			return(rez);
		},
		# fm[i, j:(j+length(value)-1)] = value
		writeSubCol = function(i, j, value) {
			stopifnot( i>=1  );
			stopifnot( j>=1  );
			stopifnot( i<=nr );
			stopifnot( j<=nc );
			stopifnot( i + length(value) - 1 <= nr );
			writeSeq( (j-1)*nr+i, value);
			return(invisible(.self));
		},
		# fm[]
		readAll = function() {
			rez = readCols(1, nc);
			return(rez);
		},
		# fm[] = value
		writeAll = function(value) {
			stopifnot( length(value) == nr*nc );
			writeSeq(1, value);
			return(invisible(.self));
		},
		appendColumns = function(mat) {
			mat = as.matrix(mat);
			oldn = .self$nc;
			newn = oldn + ncol(mat);
			if(.self$nr == 0)	{
				.self$nr = nrow(mat);
			} else {
				stopifnot( .self$nr == nrow(mat) );
			}
			.self$nc = newn;
			saveInfo();
			writeCols(oldn+1, mat);
			return(invisible(.self));
		}
	)
)

### Accessing as a usual matrix 

### Reading via vector or matrix interface.
"[.filematrix" = function(x, i, j) {
	# Basic checks
	if( !x$isOpen() )
		stop( "File matrix is not open");
	if( !missing(i) )
		if( !is.numeric(i) )
			i = round(i);
	if( !missing(j) )
		if( !is.numeric(j) )
			i = round(j);
	
	### full matrix access
	if( missing(i) & missing(j) ) {
		return( x$readAll() );
	}
	
	
	### vector access
	if(nargs()==2 & !missing(i)) {
		### checks and logical access
		if( is.logical(i) ) {
			stopifnot( length(i) == x$nr*x$nc );
			i = which(i);
		} else {
			stopifnot( min(i) >= 1 );
			stopifnot( max(i) <= x$nr*x$nc );
		}
		
		ind = .index.splitter(i);
		if(ind$n == 1) {
			rez = x$readSeq(ind$start, ind$len);
			return(rez);
		}
		rez = vector(x$type, length(i));
		for(a in 1:ind$n) {
			rez[ ind$ix[a]:(ind$ix[a+1]-1) ] =
				x$readSeq(ind$start[a], ind$len[a]);
		}
		return(rez);
	}

	### checks and logical access
	if(!missing(j)) {
		if( is.logical(j) ) {
			stopifnot( length(j) == x$nc );
			j = which(j);
		} else {
			stopifnot( min(j) >= 1 );
			stopifnot( max(j) <= x$nc );
		}
	}
	if(!missing(i)) {
		if( is.logical(i) ) {
			stopifnot( length(i) == x$nr );
			i = which(i);
		} else {
			stopifnot( min(i) >= 1 );
			stopifnot( max(i) <= x$nr );
		}
	}
	
	### column access
	if( missing(i) & !missing(j) ) {
		ind = .index.splitter(j);
		if(ind$n == 1) {
			rez = x$readCols(ind$start, ind$len);
			return(rez);
		}
		rez = vector(x$type, length(j)*x$nr);
		dim(rez) = c(x$nr, length(j));
		for(a in 1:ind$n) {
			rez[,ind$ix[a]:(ind$ix[a+1]-1)] = 
					x$readCols(ind$start[a], ind$len[a]);
		}
		return(rez);			
	}
	
	### row access via full access
	if( !missing(i) & missing(j) ) {
				j = 1:x$nc;
	}
	
	### full access
	if( !missing(i) ) {
		rez = vector(x$type, as.double(length(j)) * length(i));
		dim(rez) = c(length(i), length(j));
		if(all(diff(i)==1L)) {
			for(a in seq_along(j)) {
				rez[,a] = x$readSubCol(i[1], j[a], length(i));
			}
			return(rez);	
		} else {
			low = min(i);
			len = max(i) - low + 1;
			inew = i-low+1;
			for(a in seq_along(j)) {
				vec = x$readSubCol(low, j[a], len);
				rez[,a] = vec[inew];
			}
			return(rez);
		}
	}
	stop("What??");
}

### Writing via vector or matrix interface.
"[<-.filematrix" = function(x, i, j, value) {
	# Basic checks
	if( !x$isOpen() )
		stop( "File matrix is not open");
	if( !missing(i) )
		if( !is.numeric(i) )
			i = round(i);
	if( !missing(j) )
		if( !is.numeric(j) )
			i = round(j);
	
	### full matrix access
	if( missing(i) & missing(j) ) {
		stopifnot(length(value) == x$nr*x$nc);
		x$writeAll(value);
		return(x);
	}

	### vector access
	if(nargs()==3 & !missing(i)) {
		### checks and logical access
		if( is.logical(i) ) {
		stopifnot( length(i) == x$nr*x$nc );
		 	i = which(i);
		} else {
			stopifnot( min(i) >= 1 );
			stopifnot( max(i) <= x$nr*x$nc );
		}
		stopifnot( length(i) == length(value) );
		ind = .index.splitter(i);
		if(ind$n == 1) {
			x$writeSeq(ind$start, value);
			return(x);
		}
		for(a in 1:ind$n) {
			x$writeSeq(ind$start[a], value[ind$ix[a]:(ind$ix[a+1]-1)]);
		}
		return(x);
	}
	
	### checks and logical access
	if(!missing(j)) {
		if( is.logical(j) ) {
			stopifnot( length(j) == x$nc );
			j = which(j);
		} else {
			stopifnot( min(j) >= 1 );
			stopifnot( max(j) <= x$nc );
		}
	}
	if(!missing(i)) {
		if( is.logical(i) ) {
			stopifnot( length(i) == x$nr );
			i = which(i);
		} else {
			stopifnot( min(i) >= 1 );
			stopifnot( max(i) <= x$nr );
		}
	}	
	
	### column access
	if( missing(i) & !missing(j) ) {
		stopifnot( length(j)*x$nr == length(value) );
		dim(value) = c(x$nr, length(j));
		ind = .index.splitter(j);
		if(ind$n == 1) {
			x$writeCols(ind$start, value);
			return(x);
		}
		for(a in 1:ind$n) {
			x$writeCols(ind$start[a], value[,ind$ix[a]:(ind$ix[a+1]-1)]);
		}
		return(x);	
	}
	
	### row access via full access
	if( !missing(i) & missing(j) ) {
		j = 1:x$nc;
	}
	
	### full access
	if( !missing(i) ) {
		stopifnot( all(diff(i)==1L) );
		stopifnot( length(i)*length(j) == length(value) );
		dim(value) = c(length(i),length(j));
		for(a in seq_along(j)) {
			x$writeSubCol(i[1], j[a], value[,a]);
		}
		return(x);	
	}
	
	stop("What??");
}

### Creators of filematrix objects

# Create new, erase if exists
fm.create = function(filenamebase, nrow = 0, ncol = 1, type="double", size=NULL){
	rez = new("filematrix");
	rez$create(filenamebase=filenamebase, nrow=nrow, ncol = ncol, type=type, size=size);
	return(rez);
}

# From existing matrix
fm.create.from.matrix = function(filenamebase, mat, size=NULL) {
	rez = new("filematrix");
	rez$createFromMatrix(filenamebase=filenamebase, mat=mat, size=size);
	return(rez);
}

# Open existing file matrix
fm.open = function(filenamebase, readonly = FALSE) {
	rez = new("filematrix");
	rez$open(filenamebase=filenamebase, readonly);
	return(rez);
}

# Open and read the the whole matrix in memory.
fm.load = function(filenamebase) {
	fm = fm.open(filenamebase=filenamebase, readonly = TRUE);
	mat = as.matrix(fm);
	dimnames(mat) = dimnames(fm);
	fm$close();
	return(mat);
}

# Create from a text file matrix
fm.create.from.text.file = function(textfilename, filenamebase, skipRows = 1, skipColumns = 1, sliceSize = 1000, omitCharacters = "NA", delimiter = "\t", rowNamesColumn = 1, type="double", size = NULL) {

	s = function(x)formatC(x=x, digits=ceiling(log10(max(x)+1)), big.mark=",", big.interval=3);

	stopifnot( (skipColumns == 0) || (rowNamesColumn <= skipColumns) )
	stopifnot( (skipColumns == 0) || (rowNamesColumn >= 1) )

	fid = file(description = textfilename, open = "rt", blocking = FALSE, raw = FALSE)
	
	# clean object if file is open

	fm = fm.create(filenamebase, nrow = 1, ncol = 1, type = type, size = size);
	dim(fm) = c(0,0);
	
	lines = readLines(con = fid, n = max(skipRows,1L), ok = TRUE, warn = TRUE)
	line1 = tail(lines,1L);
	splt = strsplit(line1, split = delimiter, fixed = TRUE);
	if( skipRows > 0L ) {
		tempcolnames = splt[[1]]; # [ -(1:fileSkipColumns) ];
	} else {
		seek(fid, 0);
	}		
	
	rm( lines, line1, splt );
	
	rowNameSlices = vector("list", 15);

	curSliceId = 0L;
	repeat
	{
		# preallocate data
		curSliceId = curSliceId + 1L;
		if(length(rowNameSlices) < curSliceId) {
			rowNameSlices[[2L*curSliceId]] = NULL;
		}
		
		# read sliceSize rows
		rowtag = character(sliceSize);
		rowvals = vector("list",sliceSize);
		for(i in 1:sliceSize) {
			if( skipColumns > 0L ) {
				temp = scan(file = fid, what = character(), n = skipColumns, quiet = TRUE,sep = delimiter);
			} else {
				temp = "";
			}

			rowtag[i] = temp[rowNamesColumn];#paste(temp,collapse=" ");
			rowvals[[i]] = scan(file = fid, what = double(), nlines = 1, quiet = TRUE, sep = delimiter, na.strings = omitCharacters);
			
			if( length(rowvals[[i]]) == 0L ) {
				if(i==1L) {
					rowtag = matrix(0, 0, 0);
					rowvals = character(0);
				} else 	{
					rowtag  = rowtag[  1:(i-1) ];
					rowvals = rowvals[ 1:(i-1) ];
				}
				break;			
			}
		}
		if( length(rowtag) == 0L ) {
			curSliceId = curSliceId - 1L;
			break;
		}
		rowNameSlices[[curSliceId]] = rowtag;
		data = c(rowvals, recursive = TRUE);
		dim(data) = c(length(rowvals[[1]]), length(rowvals));
# 		data = t(data);
		fm$appendColumns(data);

		if( length(rowtag) < sliceSize ) {
			break;
		}
		cat( "Rows read: ", s(ncol(fm)), "\n");
		flush.console()
	}
	close(fid);
	if( skipRows == 0 ) {
		rownames(fm) = paste0("Col_", 1:nrow(fm));
	} else {
		rownames(fm) = tail(tempcolnames, nrow(fm));
	}
	if( skipColumns == 0 ) {
		colnames(fm) = paste0("Row_", 1:ncol(fm));
	} else {
		colnames(fm) = unlist(rowNameSlices);
	}
	cat("Rows read: ", ncol(fm), " done.\n");
	return(fm);
}


### Common interface methods

setGeneric("close")#, def = function(con){standardGeneric("close")})
setMethod("close", signature(con="filematrix"), function(con) con$close());

setGeneric("as.matrix")
setMethod("as.matrix", signature(x="filematrix"), function(x) x$readAll());

setGeneric("dim");
setMethod("dim", signature(x="filematrix"), function(x) c(x$nr, x$nc));
# dim.filematrix = function(x) as.integer(c(x$nr, x$nc));

# setGeneric("dim<-",def = function(x,value){standardGeneric("dim<-")});
setMethod("dim<-", signature(x="filematrix", value = "ANY"),	
	function(x, value) {
		x$nr = value[1];
		x$nc = value[2];
		x$saveInfo();
		return( x );
	}
);

setGeneric("length");
setMethod("length", signature(x="filematrix"),	function(x) x$nr*x$nc);

setGeneric("dimnames");
setMethod("dimnames", signature(x="filematrix"),	function(x) x$getdimnames());

setGeneric("dimnames<-");
setMethod("dimnames<-", signature(x="filematrix", value = "ANY"), function(x, value) x$setdimnames(value));

setGeneric("rownames");
setMethod("rownames", signature(x="filematrix"),	function(x) x$getrownames());

setGeneric("rownames<-");
setMethod("rownames<-", signature(x="filematrix", value = "ANY"), function(x, value) x$setrownames(value));

setGeneric("colnames");
setMethod("colnames", signature(x="filematrix"),	function(x) x$getcolnames());

setGeneric("colnames<-");
setMethod("colnames<-", signature(x="filematrix", value = "ANY"), function(x, value) x$setcolnames(value));


