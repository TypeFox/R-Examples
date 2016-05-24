#' Read a text file into a single string
#' 
#' @param filename character. Input filename.
#' @return character
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
slurp = function(filename) {
	f = file(filename, "r")
	tryCatch({
				res =  paste(readLines(f), collapse = "\n")
			}, finally = {
				close(f)
			})
	res
}

#' Write strings to a file
#' 
#' @param s character. Strings to write.
#' @param filename character. Path to output file.
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
spit = function(s, filename){
	f = file(filename, "w") 
	tryCatch({
				writeLines(con = f, text = s)
			}, finally = {
				close(f)
			})
}

#' Eval R expressions from a file.
#' 
#' @param filename character
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
evalFile = function(filename) {
	eval(parse(text = slurp(filename)))
}


#' Read big.matrix .desc file
#' 
#' @param desc_filename character. Path to .desc file
#' @return description object
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
readDesc = function(desc_filename) {
	tryCatch({
				readRDS(desc_filename)
			}, error = function(err) {
				if(err$message == "unknown input format") {
					evalFile(desc_filename)
				} else {
					stop("Bad desc file?")
				}
			})
}


#' Save big.matrix description object to disk
#' 
#' Binary format is used exclusively. 
#' 
#' @param desc_obj big.matrix description object
#' @param desc_filename character. Output file description file path.
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
saveDesc = function(desc_obj, desc_filename) {
	saveRDS(desc_obj, file = desc_filename)
#	dput(desc_obj, desc_filename)
}

#' Check equality of two lists
#' @param list1 list
#' @param list2 list
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
listEqual = function(list1, list2) {
	list1 = unlist(list1, recursive = TRUE)
	list2 = unlist(list2, recursive = TRUE)
	all(na.omit(list1 == list2))
}


type2Bytes = function(type = "double") {
	if (type == "integer")
		4
	else if (type == "double")
		8
	else if (type == "short")
		2
	else if (type == "char")
		1
	else
		stop("Unknown type!")
}


#' Coerce an R vector/matrix/data.frame into a big.matrix
#' 
#' This is a patched version of as.big.matrix from the bigmemory package.
#' The patch allows you to omit colnames/rownames even when they exist in
#' the R object.
#' 
#' 
#' @param x vector, matrix, or data.frmae
#' @param dimnames logical. FALSE by default
#' @param type See \code{bigmemory::as.big.matrix}
#' @param separated  See \code{bigmemory::as.big.matrix}
#' @param backingfile   See \code{bigmemory::as.big.matrix}
#' @param backingpath   See \code{bigmemory::as.big.matrix}
#' @param descriptorfile   See \code{bigmemory::as.big.matrix}
#' @param binarydescriptor   See \code{bigmemory::as.big.matrix}
#' @param shared   See \code{bigmemory::as.big.matrix}
#' 
#' @return big.matrix object
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @docType methods
#' @export
setGeneric('asBigMatrix', 
		function(x, type=NULL, separated=FALSE,
				backingfile=NULL, backingpath=NULL,
				descriptorfile=NULL, binarydescriptor=FALSE, shared=TRUE, 
				dimnames = FALSE) standardGeneric('asBigMatrix'))


#' @rdname asBigMatrix
#' @export 
setMethod('asBigMatrix', signature(x='matrix', dimnames = "logical"),
		function(x, type, separated, backingfile, backingpath, descriptorfile,
				binarydescriptor, shared, dimnames)
		{
			if (!is.numeric(x)) {
				warning("Casting to numeric type")
				x <- matrix(as.numeric(x), nrow=nrow(x), dimnames=dimnames(x))
			}
			if(dimnames) {
				dim_names = dimnames(x)
			} else {
				dim_names = NULL
			}
			if (is.null(type)) type <- typeof(x)
			
			if (type=="integer" | type=="double" | type=="short" | type=="char") 
			{
				y <- bigmemory::big.matrix(nrow=nrow(x), ncol=ncol(x), type=type, 
						init=NULL, dimnames=dim_names, separated=separated,
						backingfile=backingfile, backingpath=backingpath,
						descriptorfile=descriptorfile, binarydescriptor=binarydescriptor,
						shared=shared)
				y[1:nrow(x),1:ncol(x)] <- x
				junk <- gc() 
			} else stop('bigmemory: that type is not implemented.')
			return(y)
		})

#' @rdname asBigMatrix
#' @export 
setMethod('asBigMatrix', signature(x='data.frame', dimnames = "logical"),
		function(x, type, separated, backingfile, backingpath, descriptorfile,
				binarydescriptor, shared, dimnames)
		{
			warning("Coercing data.frame to matrix via factor level numberings.")
			if (is.null(type)) type <- options()$bigmemory.default.type
			if (type=="integer" | type=="double" | type=="short" | type=="char") 
			{
				if(dimnames) {
					dim_names = dimnames(x)
				} else {
					dim_names = NULL
				}
				y <- bigmemory::big.matrix(nrow=nrow(x), ncol=ncol(x), type=type, 
						init=NULL, dimnames=dim_names, separated=separated,
						backingfile=backingfile, backingpath=backingpath,
						descriptorfile=descriptorfile, binarydescriptor=binarydescriptor,
						shared=shared)
				oldbtw <- options()$bigmemory.typecast.warning
				options(bigmemory.typecast.warning=FALSE)
				for (i in 1:ncol(x)) {
					if (is.character(x[,i])) x[,i] <- factor(x[,i])
					if (is.factor(x[,i])) x[,i] <- as.numeric(x[,i])
					y[,i] <- x[,i]
				}
				options(bigmemory.typecast.warning=oldbtw)
				junk <- gc() 
			} else stop('bigmemory: that type is not implemented.')
			return(y)
		})

#' @rdname asBigMatrix
#' @export 
setMethod('asBigMatrix', signature(x='vector', dimnames = "logical"),
		function(x, type, separated, backingfile, backingpath, descriptorfile,
				binarydescriptor, shared, dimnames)
		{
			if (!is.numeric(x)) {
				warning("Casting to numeric type")
				x <- as.numeric(x)
			}
			x <- matrix(x, length(x), 1)
			warning("Coercing vector to a single-column matrix.")
			return(asBigMatrix(x, type, separated, backingfile, 
							backingpath, descriptorfile, binarydescriptor, shared, dimnames))
		})


bmBinFilename = function(stem) {
	sprintf("%s.bin", stem)
}

bmDescFilename = function(stem) {
	sprintf("%s.desc", stem)
}


#' Convert a .bin filename to a .desc filename
#' 
#' @param bin_file character. .bin filename
#' @return character
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
bin2DescFilename = function(bin_file) {
	bmDescFilename(tools::file_path_sans_ext(bin_file))
}

#' Convert a .desc filename to a .bin filename
#' 
#' @param desc_file character. .desc filename
#' @return character
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
desc2BinFilename = function(desc_file) {
	bmBinFilename(tools::file_path_sans_ext(desc_file))
}

#' Read columns into an R matrix from a big.matrix .bin file
#' 
#' @param bin_file character. Path to .bin file
#' @param ncols_to_read integer.
#' @return matrix
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
readBmBin = function(bin_file, ncols_to_read) {
	desc_file = bin2DescFilename(bin_file)
	desc = readDesc(desc_file)
	nrow = desc@description$nrow
	fh = file(bin_file, "rb")
	tryCatch({
				res = readBin(fh, 
						{
						t = desc@description$type
						if(t == "integer" || t == "double") {
							t
						} else {
							"raw"
						}
						},
						desc@description$nrow * ncols_to_read)
			}, finally = {
				close(fh)
			})
	matrix(res, nrow)
}

#' Conversion function to use when appending values to a big.matrix
#' 
#' @param desc description object
#' @return conversion function.
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
bmConvertFun = function(desc) {
	type = desc@description$type
	if(type == "double") {
		as.numeric
	} else if(type == "integer") {
		as.integer
	} else {
		stop("convert function only supports double or integer")
	}
}


#' Generate a single alpha-numeric random string
#' 
#' @param string_length integer.
#' @return character.
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
randomString = function(string_length = 6) {
	strConcat(sample(alphaNumeric, string_length))
}

#' Generate random strings
#' 
#' @param n integer. Number of string to generate.
#' @param string_length integer. Length of each string.
#' @return character.
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
randomStrings = function(n, string_length = 6) {
	sapply(1:n, function(i) randomString(string_length))
}


#' Add column(s) to an existing big.matrix
#' 
#' This function provides an effecient way to append columns to a big.matrix (
#' without copying columns that are already on disk).
#' 
#' @param bin_file character. Path to .bin file for file-backed big.matrix
#' @param dat vector, matrix or data.frame. Coercion rules are the same as in big.matrix
#' @return updated description object. 
#' @importFrom collUtils truncateEndOfFile
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
bmAddCol = function(bin_file, dat) {
	desc_file = bin2DescFilename(bin_file)
	desc = readDesc(desc_file)

	converter = bmConvertFun(desc)
	convertCol = function(vec) {
		if (is.character(vec)) vec <- factor(vec)
		if (is.factor(vec)) vec <- as.numeric(vec)
		converter(vec)
	}
	
	n_row = desc@description$nrow
	n_col = desc@description$ncol
#	print("\n==========")
#	print(class(dat))
#	print(dim(dat))
#	print(head(dat))
#	print(tail(dat))
	if(is.vector(dat)) {
		stopifnot(length(dat) == n_row)
		dat = convertCol(dat)
	} else if(is.data.frame(dat)) {
		stopifnot(nrow(dat) == n_row)
		for(i in 1:ncol(dat)) {
			dat[, i] = convertCol(dat[, i])
		}
		dat = do.call(c, dat)
	} else if(is.matrix(dat)) {
		stopifnot(nrow(dat) == n_row)
		dat = convertCol(as.vector(dat))
	}
	
	# On some systems, bin file has a trailing null byte. This is a bug in the bigmemory package.
	# I provide a temporary fix here. Bug has been reported on github.
	old_size = file.info(bin_file)$size
    n_trailing_bytes = old_size %% (n_row * n_col)
	if(n_trailing_bytes != 0) {
		collUtils::truncateEndOfFile(bin_file, n_trailing_bytes)
	}

	# write dat as new column(s)
	fh = file(bin_file, "ab")
	tryCatch({
				writeBin(dat, fh)
#				if(n_trailing_bytes != 0) {
#					writeBin(as.raw(0), fh)
#				}
			}, finally = {
				close(fh)
			})
	
	# update description
	desc = correctDesc(desc_file)
	saveDesc(desc, desc_file)

	invisible(desc)
}


#' Correct description of big.matrix 
#' 
#' @param desc_file character. Path to description file
#' @return list. Corrected description object.
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
correctDesc = function(desc_file) {
	bin_file = desc2BinFilename(desc_file)
	desc = readDesc(desc_file)
	old_n_col = desc@description$ncol
	n_row = as.numeric(desc@description$nrow)
	bytes_per_elem = type2Bytes(desc@description$type)
	bytes_per_col = n_row * bytes_per_elem
	bin_size = file.info(bin_file)$size
	new_n_col = floor(bin_size / bytes_per_col)
	if(new_n_col != old_n_col) {
		desc@description$totalCols = 
				desc@description$colOffset[2] = 
				desc@description$ncol = new_n_col
		
	}
	desc
}

#' Generate a big.matrix filename (.bin or .desc)
#' 
#' @param mat_name character. Stem of filename.
#' @param type character. Either "bin" or "desc"
#' @return character. big.matrix filename
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
bmFilename = function(mat_name, type) {
	stopifnot(type %in% c("bin", "desc"))
	sprintf("%s.%s", mat_name, type)
}

#' Get the big.matrix file path according to GCDH task tag
#' 
#' @param tag character. GCDH task tag.
#' @param mat_name character. nmiss, beta, stat, p, etc.
#' @param type character. Either "bin" or "desc"
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
bmFilepath = function(tag, mat_name, type) {
	file.path(
			tag2Dir(tag, "gcdh"),
			bmFilename(mat_name, type)
			)
}

#' Create a big.matrix under specified GCDH tag
#' @param tag character. GCDH tag.
#' @param bm_name character. Name of the big.matrix to be created.
#' @param nrow integer. Number of rows of the big.matrix
#' @param ncol integer. Number of columns of the big.matrix. Default to 1.
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
gcdhBmCreate = function(tag, bm_name, nrow, ncol = 1) {
	bmCreate(
			tag = tag, 
			type = "gcdh", 
			bm_name = bm_name, 
			nrow = nrow, 
			ncol = ncol
			)
}

bmCreate = function(tag, type, bm_name, nrow, ncol = 1) {
	backingpath = tag2Dir(tag, type)
	dir.create2(backingpath)
	bin_file = bmFilename(bm_name, "bin")
	desc_file = bmFilename(bm_name, "desc")
	res = bigmemory::filebacked.big.matrix(
			nrow = nrow, 
			ncol = ncol, 
			type = "double", 
			backingfile      = bin_file,
			descriptorfile   = desc_file,
			backingpath      = backingpath,
			binarydescriptor = TRUE
	)	
	res
}


#' Attach a big.matrix by its bin filename
#' @param bin_file character. big.matrix bin filename
#' 
#' @author Kaiyin Zhong
#' @export
bmAttachBin = function(bin_file) {
	bm = bigmemory::attach.big.matrix(bin2DescFilename(bin_file))
	bm
}



