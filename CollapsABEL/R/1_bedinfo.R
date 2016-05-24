plinkTrio <- function(bedstem, must_exist = FALSE) {
	ext_trio = c("bed", "bim", "fam")
	plink_trio = paste(bedstem, ext_trio, sep = ".")
	if(must_exist) {
		filePath(plink_trio)@path
	} else {
		plink_trio
	}
}

#' An S4 class representing info about plink files
#' 
#' Info about plink files, including the root directory, 
#' paths of plink .bed, .bim, .fam and .frq files, ff backing 
#' directories for .bim, .fam and .frq files, etc.
#' 
#' @slot main_dir Root directory where .bed, .bim and .fam files sit.
#' @slot plink_stem character. Path to the .bed file sans the extension name
#' @slot plink_trio character of length 3. Paths to .bed, .bim and .fam files (in that order).
#' @slot plink_trio_base character. Basenames of \code{plink_trio}.
#' @slot plink_frq character. Path to .frq file.
#' 
#' @export 
.PlInfoC = setClass("PlInfoC", representation(
				main_dir = "character", 
				plink_stem = "character",
				plink_trio = "character", 
				plink_trio_base = "character", 
				plink_frq = "character"
		), 
		prototype(
				main_dir = "", 
				plink_stem = "",
				plink_trio = rep("", 3), 
				plink_trio_base = rep("", 3), 
				plink_frq = ""
		), validity = function(object) {
			obj_slots = list(
					object@main_dir,
					object@plink_stem ,
					object@plink_trio ,
					object@plink_trio_base ,
					object@plink_frq 
			) 
			names(obj_slots) = c(
					"main_dir",
					"plink_stem",
					"plink_trio",
					"plink_trio_base",
					"plink_frq"
			)
			msg = lenCheck(
					obj_slots,
					c(1, 1, 3, 3, 1))			
			if(msg != TRUE) {
				return(msg)
			}
			ext_trio = c("bed", "bim", "fam")
			if(!all(tools::file_ext(object@plink_trio) 
							== ext_trio)) {
				return(paste("Extensions should be: ", strConcat(ext_trio)))
			}
			miss_files = nonExistentFiles(object@plink_trio)
			if(length(miss_files) > 0) {
				return(miss_files)
			} else {
				return(TRUE)
			}
		})


#' Constructor for PlInfoC class
#' 
#' Populates an PlInfoC object from a given plink bed filename stem (i.e. exclude extension name)
#' 
#' @param pl_info a PlInfoC object, possibly empty.
#' @param bedstem path of bed file excluding extension name
#' @param db_setup logical. Whether to setup SQLite database for .bim, .fam and .frq files.
#' @return a PlInfoC object
#' @examples 
#' \dontrun{
#'			pl_info = plInfo(.PlInfoC(), "mmp13", db_setup = TRUE)
#'			isSetup(pl_info)
#'			bim_ff = getQuery(sqliteFilePl(pl_info), "select * from bim")
#'			fam_ff = getQuery(sqliteFilePl(pl_info), "select * from fam")
#'			frq_ff = getQuery(sqliteFilePl(pl_info), "select * from frq")
#' }
#' @importFrom methods validObject
#' @author Kaiyin Zhong, Fan Liu
#' @export 
setGeneric("plInfo",
		function(pl_info, bedstem, db_setup) {
			standardGeneric("plInfo")
		})

#' @rdname plInfo
#' @export 
setMethod("plInfo",
		signature(pl_info = "PlInfoC", bedstem = "character", db_setup = "logical"),
		function(pl_info, bedstem, db_setup) { 
			# plink trio
			ext_trio = c("bed", "bim", "fam")
			plink_trio = normalizePath(
					plinkTrio(bedstem = bedstem, must_exist = FALSE)
			) 
			plink_trio_base = basename(plink_trio)
			names(plink_trio) = names(plink_trio_base) = ext_trio
			plink_stem = tools::file_path_sans_ext(plink_trio["bed"])
			names(plink_stem) = NULL
			
			# main dir where plink files sit
			main_dir = dirname(plink_trio[1])
			
			# frq file
			plink_frq = paste(bedstem, ".frq", sep="")
			
			# return a PlInfoC obj
			pl_info@main_dir = main_dir
			pl_info@plink_stem = plink_stem
			pl_info@plink_trio = plink_trio
			pl_info@plink_trio_base = plink_trio_base
			pl_info@plink_frq = plink_frq
			validObject(pl_info)
			if(db_setup) {
				setup(pl_info)
			}
			pl_info
		})

#' @rdname plInfo
#' @export 
setMethod("plInfo",
		signature(pl_info = "PlInfoC", bedstem = "character", db_setup = "missing"),
		function(pl_info, bedstem, db_setup) { 
			plInfo(pl_info, bedstem, FALSE)
		})

#' @rdname plInfo
#' @export 
setMethod("plInfo",
		signature(pl_info = "missing", bedstem = "character", db_setup = "logical"),
		function(pl_info, bedstem, db_setup) {
			plInfo(.PlInfoC(), bedstem, db_setup)
		})

#' @rdname plInfo
#' @export 
setMethod("plInfo",
		signature(pl_info = "missing", bedstem = "character", db_setup = "missing"),
		function(pl_info, bedstem, db_setup) {
			plInfo(.PlInfoC(), bedstem, FALSE)
		})

#' SQLite file of a PlInfoC object
#' 
#' @param x PlInfoC or PlGwasC object
#' @return character. Path to SQLite database file.
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
sqliteFilePl = function(x) {
	if(isS4Class(x, "PlInfoC")) {
		filename = sprintf("%s.sqlite", x@plink_stem)
	} else if(isS4Class(x, "PlGwasC")) {
		filename = sprintf("%s.sqlite", x@pl_info@plink_stem)
	}
	filename
}

#' Check if a directory containing .bed .fam and .bim files is properly setup
#'  
#' @param pl_info PlInfoC object
#' @return TRUE or FALSE
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
isSetup = function(pl_info) {
	stopifnot(isS4Class(pl_info, "PlInfoC"))
	sql_file = sqliteFilePl(pl_info)
	isSQLite3(sql_file) && dbUpToDate(sql_file)
}

#' Setup up a directory containing plink files 
#' 
#' @param pl_info PlInfoC object
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
setup = function(pl_info) {
	stopifnot(isS4Class(pl_info, "PlInfoC"))
	if(isSetup(pl_info)) {
		TRUE
	} else {
		sqlite_file = sqliteFilePl(pl_info)
		if(file.exists(sqlite_file)) {
			file.remove(sqlite_file)
		}
		if(!file.exists(pl_info@plink_frq) || !frqUpToDate(pl_info@plink_frq)) {
			plinkr(bfile = pl_info@plink_stem, 
					freq = "", 
					out = pl_info@plink_stem, 
					wait = TRUE)
		}
		frq = read.table(pl_info@plink_frq, header = TRUE, stringsAsFactors = FALSE)
		bim = readBim(pl_info@plink_trio["bim"])
		fam = readFam(pl_info@plink_trio["fam"])
		fam = setNames(fam, c("FID", "IID", "PID", "MID", "SEX", "PHE"))
		tryCatch({
					file.create2(sqlite_file)
					db = RSQLite::dbConnect(RSQLite::SQLite(), sqlite_file)
					RSQLite::dbWriteTable(db, "bim", bim)
					RSQLite::dbWriteTable(db, "fam", fam)
					RSQLite::dbWriteTable(db, "frq", frq)
				}, finally = {
					RSQLite::dbDisconnect(db)
				})			
	}
}


#' Get number of individuals
#' 
#' @param pl_info PlInfoC object
#' @export 
nIndivPl = function(pl_info) {
	stopifnot(isS4Class(pl_info, "PlInfoC"))
	getQuery(sqliteFilePl(pl_info), "select count(iid) from fam")[1, 1]
}

#' Get number of SNPs.
#' 
#' @param pl_info PlInfoC object
#' @export 
nSnpPl = function(pl_info) {
	stopifnot(isS4Class(pl_info, "PlInfoC"))
	getQuery(sqliteFilePl(pl_info), "select count(snp) from bim")[1, 1]
}

#' Get number of bytes used by each SNP.
#' 
#' @param pl_info PlInfoC object
#' @export 
bytesSnp = function(pl_info) {
	stopifnot(isS4Class(pl_info, "PlInfoC"))
	as.numeric(ceiling(nIndivPl(pl_info) / 4))
}

#' Get apparent number of individuals
#' 
#' @param pl_info PlInfoC object
#' @export 
nIndivApprPl = function(pl_info) {
	stopifnot(isS4Class(pl_info, "PlInfoC"))
	as.numeric(bytesSnp(pl_info) * 4)
}


#' FID and IID columns from fam file
#' 
#' @param pl_info PlInfoC object
#' @return data.frame of two columns "FID" and "IID"
#' 
#' @examples 
#' \dontrun{
#' pl_info = plInfo(bedstem = "mmp13", db_setup = TRUE)
#' fidiid = fidIid(pl_info)
#' fam = readFam("mmp13.fam", c("FID", "IID"))
#' all(fam == fidiid)
#' }
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
fidIid = function(pl_info) {
	setup(pl_info)
	getQuery(sqliteFilePl(pl_info), "select fid, iid from fam order by rowid")
}
