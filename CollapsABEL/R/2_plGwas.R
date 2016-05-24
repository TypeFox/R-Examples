#' An S4 class representing info about GWAS on plink files
#' 
#' @slot gwas_tag character. Tag for this GWAS.
#' @slot opts list. Plink options.
#' @importFrom stringr str_match 
#' @importFrom dplyr left_join
#' 
#' @export 
.PlGwasC = setClass("PlGwasC", 
		representation(
				gwas_tag = "character",
				opts = "list"
		), 
		prototype(
				gwas_tag = "",
				opts = list()
		), 
		contains = "RbedInfoC", 
		validity = function(object) {
			obj_slots = list( object@gwas_tag)
			names(obj_slots) = c( "gwas_tag")
			msg = lenCheck(obj_slots, 1)
			if(msg != TRUE) {
				return(msg)
			}
			
			if(!file.exists(object@opts$pheno)) {
				return("plink phenotype file does not exist.")
			}
			
			first_line = readLiteral(object@opts$pheno, nrows = 1)[1, ]
			
			if(length(object@opts$pheno_name) != 1) 
				return("one phenotype at a time.")
			
			if(!(object@opts$pheno_name %in% first_line)) {
				return(
						sprintf("phenotype name (%s) not in phenotype file.", 
								object@opts$pheno_name))
			}
			
			covars = covarNames(object)
			if(covars[1] != "") {
				covar_not_there = which(! covars %in% first_line)
				if(length(covar_not_there) > 0) {
					return(
							strConcat(
									c(
											"Covars not found: ",
											strConcat(covars[covar_not_there], " ")
									)))
				}
			}
			
			TRUE
		})



#' Get covariate names of a GWAS
#' 
#' @name covarNames
#' 
#' @param pl_gwas PlGwasC object.
#' @return character. Vector of covariate names.
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @docType methods
#' @export
covarNames = function(pl_gwas) {
	stopifnot(isS4Class(pl_gwas, "PlGwasC"))
	if("covar_name" %in% names(pl_gwas@opts)) 
		strsplit(pl_gwas@opts$covar_name, ",")[[1]]
	else 
		""
}

#' Constructor for PlGwasC class
#' 
#' @param pl_gwas PlGwasC or PlInfoC object
#' @param pheno character. Phenotype file
#' @param pheno_name character. Phenotype names.
#' @param covar_name character. Covariate names.
#' @param gwas_tag character. Tag for this GWAS.
#' @param assoc logical. Whether use the "--assoc" option for PLINK.
#' @param opts list. Options to be passed to PLINK.
#' @return PlGwasC object 
#' @examples 
#' \dontrun{
#' gwas_tag = "mmp13_page_sex_age"
#' rbed_info = rbedInfo(bedstem = "mmp13")
#' pl_gwas = plGwas(rbed_info, 
#' 		pheno = "mmp13.phe",
#' 		pheno_name = "Page", 
#' 		covar_name = "Sex,Cage", 
#' 		gwas_tag = gwas_tag)
#' runGwas(pl_gwas)
#' "mmp13_page_sex_age" %in% listGwasTags() == "TRUE"
#' gwas_out = readGwasOut(pl_gwas, rmGwasOut = FALSE)
#' colClasses(gwas_out) == c("integer", "character", "integer", 
#'     "character", "character", "integer", 
#'     "numeric", "numeric", "numeric")
#' }
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @docType methods
#' @export
setGeneric("plGwas",
		function(pl_gwas, pheno, pheno_name, covar_name, gwas_tag, 
				assoc, opts) {
			standardGeneric("plGwas")
		})

#' @rdname plGwas
#' @export 
setMethod("plGwas",
		signature(pl_gwas = "PlGwasC", pheno = "character", 
				pheno_name = "character", covar_name = "character", 
				gwas_tag = "character",  
				assoc = "logical", opts = "list"
		),
		function(pl_gwas, pheno, 
				pheno_name, covar_name, 
				gwas_tag, 
				assoc, opts
		) {
			if(gwas_tag %in% listTags("gwas")) {
				removeTag(gwas_tag, "gwas")
			}
			pl_gwas@gwas_tag = gwas_tag
			pl_gwas@opts$bfile = pl_gwas@pl_info@plink_stem
			pl_gwas@opts$pheno = pheno
			pl_gwas@opts$pheno_name = pheno_name
			if(covar_name != "") {
				pl_gwas@opts$covar = pheno
				pl_gwas@opts$covar_name = covar_name			
			}
			pl_gwas@opts$allow_no_sex = ""
			pl_gwas@opts$out = gwasOutStem(pl_gwas)
			pl_gwas@opts$wait = TRUE
			if(assoc) {
				pl_gwas@opts$assoc = ""
			} else if(binPhe(pl_gwas)) {
				pl_gwas@opts$logistic = "hide-covar beta"
			} else {
				pl_gwas@opts$linear = "hide-covar"
			}
			pl_gwas
		})

#' @rdname plGwas
#' @export 
setMethod("plGwas",
		signature(pl_gwas = "RbedInfoC", pheno = "character", 
				pheno_name = "character", covar_name = "character", 
				gwas_tag = "character", 
				assoc = "logical", opts = "list"),
		function(pl_gwas, 
				pheno, pheno_name, covar_name, 
				gwas_tag, 
				assoc, opts) {
			pl_gwas = as(pl_gwas, "PlGwasC")
			plGwas(pl_gwas,  
					pheno, pheno_name, covar_name, 
					gwas_tag ,
					assoc, opts
			)
		}
)


#' @rdname plGwas
#' @export 
setMethod("plGwas",
		signature(pl_gwas = "RbedInfoC", pheno = "character", 
				pheno_name = "character", covar_name = "character", 
				gwas_tag = "character",
				assoc = "missing", opts = "missing"),
		function(pl_gwas, 
				pheno, pheno_name, covar_name, 
				gwas_tag, 
				assoc, opts) {
			plGwas(pl_gwas,  
					pheno, pheno_name, covar_name, 
					gwas_tag, 
					FALSE, list())
		}
)

#' @rdname plGwas
#' @export 
setMethod("plGwas",
		signature(pl_gwas = "RbedInfoC", pheno = "character", 
				pheno_name = "character", covar_name = "character", 
				gwas_tag = "character", 
				assoc = "missing", opts = "list"
		),
		function(pl_gwas, 
				pheno, pheno_name, covar_name, 
				gwas_tag, 
				assoc, opts) {
			plGwas(pl_gwas,  
					pheno, pheno_name, covar_name, 
					gwas_tag, 
					FALSE, opts)
		}
)

#' @rdname plGwas
#' @export 
setMethod("plGwas",
		signature(pl_gwas = "RbedInfoC", pheno = "character", 
				pheno_name = "character", covar_name = "character", 
				gwas_tag = "character", 
				assoc = "logical", opts = "missing"
		),
		function(pl_gwas, 
				pheno, pheno_name, covar_name, 
				gwas_tag, 
				assoc, opts) {
			plGwas(pl_gwas,  
					pheno, pheno_name, covar_name, 
					gwas_tag, 
					assoc, list())
		}
)

#' @rdname plGwas
#' @export 
setMethod("plGwas",
		signature(pl_gwas = "RbedInfoC", pheno = "character", 
				pheno_name = "character", covar_name = "missing", 
				gwas_tag = "character", 
				assoc = "missing", opts = "missing"
		),
		function(pl_gwas, 
				pheno, pheno_name, covar_name, 
				gwas_tag, 
				assoc, opts) {
			plGwas(pl_gwas,  
					pheno, pheno_name, "", 
					gwas_tag, 
					FALSE, list())
		}
)

setGeneric("gwasDir",
		function(pl_gwas, ...) {
			standardGeneric("gwasDir")
		})

#' GWAS results directory of a certain GWAS scan
#'  
#' @param pl_gwas PlGwasC object
#' @return character. 
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
gwasDir = function(pl_gwas) {
	stopifnot(isS4Class(pl_gwas, "PlGwasC"))
	gwas_dir = file.path(collenv$.collapsabel_gwas, pl_gwas@gwas_tag)
	if(!file.exists(gwas_dir)) {
		dir.create(gwas_dir, recursive = TRUE)
	}
	gwas_dir
}



#' Plink output filename
#' 
#' To be passed as the \code{--out} option to plink.
#' 
#' @name gwasOutStem
#' 
#' @param pl_gwas PlGwasC object.
#' @return character. Plink output filename, without extension
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
gwasOutStem = function(pl_gwas) {
	stopifnot(isS4Class(pl_gwas, "PlGwasC"))
	file.path(gwasDir(pl_gwas), basename(pl_gwas@pl_info@plink_stem))
}

#' Check whether an S4 object is of a certain class
#' 
#' @param obj S4 object
#' @param c Class name
#' @return  logical
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
isS4Class = function(obj, c) {
	isS4(obj) && is(obj, c)
}

#' GWAS output file name
#' 
#' @param pl_gwas PlGwasC object.
#' @return character
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
gwasOut = function(pl_gwas) {
	stem = gwasOutStem(pl_gwas)
	if("assoc" %in% names(pl_gwas@opts)) {
		if(binPhe(pl_gwas)) {
			paste(stem, ".assoc", sep = '')
		} else {
			paste(stem, ".qassoc", sep = '')
		}
	} else if("linear" %in% names(pl_gwas@opts)) {
		paste(stem, ".assoc.linear", sep = "")
	} else if("logistic" %in% names(pl_gwas@opts)) {
		paste(stem, ".assoc.logistic", sep = "")
	} else {
		stop("no modeling option specified?")
	}
}


#' Set analysis model
#' 
#' @param pl_gwas PlGwasC object.
#' @param mod character. One of "linear", "logistic" or "assoc", default to "linear". 
#' @return PlGwasC object
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
setOptModel = function(pl_gwas, mod = "linear") {
	poss_mods = c("linear", "logistic", "assoc")
	stopifnot(isS4Class(pl_gwas, "PlGwasC"))
	stopifnot((length(mod) == 1) && (mod %in% poss_mods))
	# nullify all modes that are not selected
	mods_to_null = poss_mods[poss_mods != mod]
	for(m in mods_to_null) {
		pl_gwas@opts[[m]] = NULL
	}
	if(mod == "logistic" || mod == "linear") {
		pl_gwas@opts[[mod]] = "hide-covar"
	} else {
		pl_gwas@opts[[mod]] = ""
		pl_gwas@opts$covar_name = NULL
	}
	pl_gwas
}


#' Run a GWAS
#' 
#' @param pl_gwas PlGwasC object
#' @param wait logical. Wait until GWAS is finished if this is set to TRUE. Default to FALSE.
#' @param save_pl_gwas logical. Whether to save the plGwas object. Default to FALSE.
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
runGwas = function(pl_gwas, wait = TRUE, save_pl_gwas = FALSE) {
	stopifnot(isS4Class(pl_gwas, "PlGwasC"))
	gwas_dir = gwasDir(pl_gwas)
	if(file.exists(gwas_dir)) {
		unlink(gwas_dir, recursive = TRUE)
	}
	dir.create2(gwas_dir)
	if(!wait) {
		pl_gwas@opts$wait = FALSE
	}
	do.call(plinkr, pl_gwas@opts)
	if(save_pl_gwas) {
		saveRDS(pl_gwas, gwasRDS(pl_gwas))
	}
}

#' Read GWAS output from plink

#' If the GWAS is finished, returns a data.frame, 
#' otherwise returns NULL.
#' @param pl_gwas PlGwasC object.
#' @param cn_select Colnames to select. Default to "..all"
#' @param rmGwasOut Logical. Whether to remove GWAS output files after finished reading them. Default to TRUE. 
#' @return data.frame 
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
readGwasOut = function(pl_gwas, cn_select = "..all", rmGwasOut = TRUE) {
	filepath = gwasOut(pl_gwas)
	res = readPlinkOut(filepath, cn_select)
	if(rmGwasOut) file.remove(filepath)
	res
}

#' Trim plink files
#' 
#' This function calculates number of individuals in .fam file (n1)
#' and number of individuals in phenotype file (n2). If n1 > n2, then
#' all the individuals not included in the phenotype file will be 
#' removed from plink files.
#' 
#' @param pl_gwas PlGwasC object.
#' @param suffix character. Suffix to the new plink file names.
#' @return PlGwasC object
#' @importFrom collUtils countlines
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
plTrim = function(pl_gwas, suffix="trimmed") {
	old_stem = pl_gwas@pl_info@plink_stem
	new_stem = paste(old_stem, suffix, sep = "_")
	if(pl_gwas@nindiv > countlines(pl_gwas@opts$pheno)) {
		plinkr(bfile = old_stem, 
				keep = pl_gwas@opts$pheno, 
				out = new_stem, 
				wait = TRUE, 
				make_bed = ""
				)
		pl_gwas_trimmed = plGwas(rbedInfo(new_stem), 
				pl_gwas@opts$pheno, 
				pl_gwas@opts$pheno_name, 
				pl_gwas@opts$covar_name, 
				paste(pl_gwas@gwas_tag, suffix, sep = "_"), 
				ifelse("assoc" %in% names(pl_gwas@opts), TRUE, FALSE)
		)
		pl_gwas_trimmed@opts$bfile = new_stem
		pl_gwas_trimmed
	} else {
		pl_gwas
	}
}

#' Remove GWAS results by tag
#' 
#' @param x character. Tag name.
#' @param type character. Type of tag.
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
removeTag = function(x, type = "gwas") {
	unlink(tag2Dir(x, type), recursive = TRUE, 
			force = TRUE)
}

removeGwasTag = removeTag

removeAllTags = function(type = "gwas") {
	tags = listTags(type = type) 
	for(tag in tags) {
		removeTag(tag, type = type)
	}
}


#' Get RDS file path of a PlGwasC object
#' 
#' @param pl_gwas PlGwasC object.
#' @return character. path of a PlGwasC object
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
gwasRDS = function(pl_gwas) {
	tag2RDSPath(pl_gwas@gwas_tag)
}


#' List GWAS or GCDH tags
#' 
#' @param type character. Either "gwas" or "gcdh".
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
listGwasTags = function(type = "gwas") {
	list.files(file.path(collenv$.collapsabel_dir, type))
}


#' @rdname listGwasTags
#' @export
listTags = listGwasTags


tag2Dir = function(x, type = "gwas") {
	file.path(collenv$.collapsabel_dir, type, x)
}


tag2RDS = function(x, type = "gwas") {
	if(type == "gwas") {
		paste(x, 
				"_PlGwasC.rds", 
				sep = "")
	} else {
		paste(x, 
				"_GCDH.rds", 
				sep = "")
	}
}

tag2RDSPath = function(x, type = "gwas") {
	d = tag2Dir(x, type)
	file.path(d, tag2RDS(x, type))
}

dir2Tag = function(x) {
	basename(x)
}


#' Load PlGwasC object by tag, from the RDS file
#' 
#' @param gwas_tag character. Tag of a GWAS run.
#' @return PlGwasC object.
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
loadGwas = function(gwas_tag) {
	stopifnot(is.character(gwas_tag))
	stopifnot(length(gwas_tag) == 1)
	gwas_dirs = listGwasTags()
	if(gwas_tag %in% gwas_dirs) {
		rds_file = tag2RDSPath(gwas_tag)
		if(file.exists(rds_file)) {
			readRDS(rds_file)
		} else {
			NULL
		}
	} else {
		NULL
	}
}



#' Get file size
#' 
#' @param filename character. Path to file.
#' @return integer. Size of file.
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
fileSize = function(filename) {
	stopifnot(is.character(filename))
	stopifnot(length(filename) == 1)
	filename = filePath(filename)@path
	file.info(filename)$size
}




#' Read phenotype file
#' 
#' @param pl_gwas PlGwasC object
#' @param cn_select Colnames to select. Default to "..all", which means all columns are read in.
#' @return data.frame
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
readPhe = function(pl_gwas, cn_select = "..all") {
	stopifnot(isS4Class(pl_gwas, "PlGwasC"))
	info = readInfo(pl_gwas@opts$pheno)
	phe = info@read_fun(info, cn_select)
	if(length(cn_select) == 1 && cn_select == "..all") {
		classes = colClasses(headPhe(pl_gwas, 1))
	} else {
		classes = colClasses(headPhe(pl_gwas, 1)[, cn_select, drop = FALSE])
	}
	phe = correctTypes(phe, types = classes)
	phe = charify(phe, c("FID", "IID"))
	phe
}

#' Convert certain columns of a data.frame to character type
#' 
#' @param dat data.frame
#' @param cols character. Names of columns to be converted.
#' @return  data.frame
#' @examples 
#' \dontrun{
#' x = data.frame(x = 1:3, y= 2:4)
#' all(colClasses(x) == c("integer", "integer"))
#' x = charify(x, "x")
#' all(colClasses(x) == c("character", "integer"))
#' }
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
charify = function(dat, cols) {
	cnames = colnames(dat)
	for(cols_i in cols) {
		if(cols_i %in% cnames) {
			dat[, cols_i] = as.character(dat[, cols_i])
		}
	}
	dat
}

#' Read first n lines of a phenotype file
#' 
#' @param pl_gwas PlGwasC object
#' @param nrows number of lines to read
#' @return data.frame
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
headPhe = function(pl_gwas, nrows = 5L) {
	stopifnot(isS4Class(pl_gwas, "PlGwasC"))
	stopifnot(is.numeric(nrows) && length(nrows) == 1)
	read.table(pl_gwas@opts$pheno, header = TRUE, nrows = nrows, 
			stringsAsFactors = FALSE)
}

#' Check whether a trait is binary
#' 
#' @param v numeric vector.
#' @param na_value a vector of numeric values which should be seen as NA.
#' @return logical
#' @examples
#' \dontrun{
#' !isBinary(c(1, 1.1, 1, 1.1, NA))
#' isBinary(c(1, 2, 1, 2, NA))
#' !isBinary(c(-9, 2.3, 4.1, -9, -9), -9)
#' isBinary(c(-9, 2, 4, -9, -9), -9)
#' isBinary(c(1, 2, 2, 1, -9, -9.9), c(-9, -9.9))
#' }
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
isBinary = function(v, na_value = NULL) {
	if(! is.numeric(v)) {
		return(FALSE)
	}
	v = unique(na.omit(v))
	if(! is.null(na_value)) {
		v = v[! v %in% na_value]
	}
	if(length(v) > 2L || any(as.integer(v) != v)) {
		FALSE
	} else {
		TRUE
	}
}

#' Check whether phenotype of a GWAS is binary
#' 
#' @param pl_gwas PlGwasC object.
#' @return logical
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
binPhe = function(pl_gwas) {
	stopifnot(isS4Class(pl_gwas, "PlGwasC"))
	validObject(pl_gwas)
	phe = readPhe(pl_gwas, pl_gwas@opts$pheno_name)
	isBinary(phe[, 1])
}

#' Plink log fie
#' 
#' Redirect stdout to this file when plink is running.
#' 
#' @param pl_gwas PlGwasC object.
#' @return character. Path to log file.
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
gwasLog = function(pl_gwas) {
	stopifnot(isS4Class(pl_gwas, "PlGwasC"))
	sprintf("%s.stdout_log", gwasOutStem(pl_gwas))
}


#' Invoke a GWAS in R
#' 
#' @param pl_gwas PlGwasC object.
#' @param snp_vec numeric or character. Vector of SNPs.
#' @return matrix. Coefficient matrix. One row for each SNP.
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
gwasR = function(pl_gwas, snp_vec) {
	pheno_name = pl_gwas@opts$pheno_name
	covars = covarNames(pl_gwas)
	gwas_dat = gwasDat(pl_gwas, snp_vec)
	if(binPhe(pl_gwas)) {
		model_family = binomial
	} else {
		model_family = gaussian
	}
	glmIter(gwas_dat$dat, 
			y = pheno_name, 
			xs = gwas_dat$geno_names, 
			covars = covars 
			)
}


#' Read genotype and phenotype data into R
#' 
#' @param pl_gwas PlGwasC object.
#' @param snp_vec numeric or character. Vector of SNPs.
#' @return data.frame
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
gwasDat = function(pl_gwas, snp_vec) {
	stopifnot(isS4Class(pl_gwas, "PlGwasC"))
	stopifnot(is.numeric(snp_vec) || is.character(snp_vec))
	geno = readBed(pl_gwas, snp_vec)
	phe = readPhe(pl_gwas)
	list(dat = dplyr::left_join(phe, geno, by = c("FID", "IID")), 
			geno_names = colnames(geno)[3:ncol(geno)],
			pheno_names = colnames(phe)[3:ncol(phe)]
	)
}



chGwasTag = function(pl_gwas, new_tag) {
	removeTag(pl_gwas@gwas_tag, "gwas")
	pl_gwas@gwas_tag = new_tag
	pl_gwas@opts$out = gwasOutStem(pl_gwas)
	pl_gwas
}


