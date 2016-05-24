#' S4 class for necessary info to read a bed file into R
#'
#' @slot pl_info PlInfoC object
#' @slot jbed jobjRef object, of Bed class in java
#' @slot nsnp numeric. Number of SNPs.
#' @slot nindiv numeric. Number of individuals.
#' @slot nindiv_appr numeric. Apparent number of individuals.
#' @slot bytes_snp numeric. Number of bytes used for each SNP.
#' @export
.RbedInfoC = setClass("RbedInfoC",
		representation(
				pl_info = "PlInfoC",
				jbed = "jobjRef",
				nsnp = "numeric",
				nindiv = "numeric",
				nindiv_appr = "numeric",
				bytes_snp = "numeric"
		),
		prototype(
				pl_info = .PlInfoC(),
				jbed = .jnew("java/lang/Integer", 0L),
				nsnp = 0,
				nindiv = 0,
				nindiv_appr = 0,
				bytes_snp = 0
		),
		validity = function(object) {
			if(! .jinstanceof(object@jbed, "vu/co/kaiyin/Bed")) {
				return("jbed is not of java Bed class.")
			}
			TRUE
		})

#' Constructor of RbedInfoC class
#'
#' @param bedstem character. Path to bed file without extension.
#' @param db_setup logical. Whether to setup SQLite database for .bim, .fam and .frq files.
#' @return An RbedInfoC object.
#' @importFrom collUtils rBed
#'
#' @author Kaiyin Zhong, Fan Liu
#' @export
rbedInfo = function(bedstem, db_setup = FALSE) {
	stopifnot(length(bedstem) == 1)
	stopifnot(is.character(bedstem))
	pl_info = plInfo(bedstem = bedstem, db_setup = db_setup)
	bed_path = pl_info@plink_trio["bed"]
	
	rbed_info = .RbedInfoC()
	rbed_info@pl_info = pl_info
	rbed_info@jbed = collUtils::rBed(bed_path)
	if(db_setup) {
		setupRbed(rbed_info)
	} else {
		rbed_info
	}
}

#' Check if an RbedInfoC object is properly set up
#' @param rbed_info RbedInfoC object
#' @return logical.
#'
#' @author Kaiyin Zhong, Fan Liu
#' @export
isSetupRbed = function(rbed_info) {
	res = isSetup(rbed_info@pl_info) &&
			rbed_info@bytes_snp > 0 &&
			rbed_info@nindiv_appr > 0 &&
			rbed_info@nindiv > 0 &&
			rbed_info@nsnp > 0
	if(is.na(res) || is.null(res)) {
		res = FALSE
	}
	res
}

#' Setup an RbedInfoC object
#'
#' The setup job includes the following tasks: 1. Set up the PlInfoC object.
#' 2. Calculate number of bytes used by each SNP. 3. Calculate the Number of individuals.
#' 4. Calculate total number of SNPs. 5. Validate the RbedInfoC object.
#'
#' @param rbed_info RbedInfoC object
#' @return RbedInfoC object
#'
#' @author Kaiyin Zhong, Fan Liu
#' @export
setupRbed = function(rbed_info) {
	if(isSetupRbed(rbed_info)) {
		rbed_info
	} else {
		pl_info = rbed_info@pl_info
		setup(pl_info)
		rbed_info@bytes_snp = bytesSnp(pl_info)
		rbed_info@nindiv = nIndivPl(pl_info)
		rbed_info@nindiv_appr = nIndivApprPl(pl_info)
		rbed_info@nsnp = nSnpPl(pl_info)
		validObject(rbed_info)
		rbed_info
	}
}


#' Check whether bed file is of correct size
#'
#' It is correct if its real size is the equal to its theoretical size.
#'
#' @name bedSizeCorrect
#'
#' @param rbed_info RbedInfoC object
#' @return logical.
#'
#' @author Kaiyin Zhong, Fan Liu
#' @export
bedSizeCorrect = function(rbed_info) {
	stopifnot(isS4Class(rbed_info, "RbedInfoC"))
	realBedSize(rbed_info) == theoBedSize(rbed_info)
}


#' File size of bed file
#'
#' @name realBedSize
#'
#' @param rbed_info RbedInfoC object
#' @return numeric. Size of bed file.
#'
#' @author Kaiyin Zhong, Fan Liu
#' @export
realBedSize = function(rbed_info) {
	stopifnot(isS4Class(rbed_info, "RbedInfoC"))
	as.numeric(file.info(rbed_info@pl_info@plink_trio["bed"])$size)
}



#' Theoretical size of bed file
#'
#' Computed from dimensions of bim an fam files.
#'
#' @name theoBedSize
#'
#' @param rbed_info RbedInfoC object
#' @return numeric. Theoretical size of bed file.
#'
#' @author Kaiyin Zhong, Fan Liu
#' @export
theoBedSize = function(rbed_info) {
	stopifnot(isS4Class(rbed_info, "RbedInfoC"))
	(as.numeric(rbed_info@nindiv_appr) *
				as.numeric(rbed_info@nsnp) / 4) + 3
}


#' Read genotypes from PLINK bed file into R
#'
#' @name readBed
#'
#' @param rbed_info RbedInfoC object
#' @param snp_vec numeric. Vector of SNP index. Either row numbers in the bim file or a vector of SNP names.
#' @param fid_iid logical. Whether the FID and IID columns should be included.
#' @param snp_names_as_colnames logical. Whether SNP names should be used as colnames in the returned data frame
#' @return data.frame Genotype data from bed file.
#' @importFrom collUtils getJArray
#'
#' @author Kaiyin Zhong, Fan Liu
#' @docType methods
#' @export
setGeneric("readBed",
		function(rbed_info, snp_vec,
				fid_iid = TRUE, snp_names_as_colnames = TRUE) {
			standardGeneric("readBed")
		})

#' @rdname readBed
#' @export
setMethod("readBed",
		signature(rbed_info = "RbedInfoC", snp_vec = "ANY",
				fid_iid = "logical", snp_names_as_colnames = "logical"),
		function(rbed_info, snp_vec,
				fid_iid, snp_names_as_colnames
		) {
#			if(is.numeric(snp_vec)) {
#				snp_vec = sort(as.integer(snp_vec))
#				snp_names = getQuery(
#						sqliteFilePl(rbed_info@pl_info),
#						sprintf("select snp from bim where rowid in %s order by rowid",
#								numVectorSQLRepr(snp_vec)))[, 1]
#			} else if(is.character(snp_vec)) {
#				snp_names = snp_vec
#				snp_vec = snpRowId(rbed_info@pl_info, snp_vec)
#				snp_ord = order(snp_vec)
#				snp_vec = snp_vec[snp_ord]
#				snp_names = snp_names[snp_ord]
			if(is.numeric(snp_vec)) {
				snp_vec = sort(as.integer(snp_vec))
			} else if(is.character(snp_vec)) {

			} else {
				stop("snp_vec must be either numeric or character.")
			}
			
#			mat_ref = rbed_info@jbed$readBed( .jarray(snp_vec))
			geno_data_ref = rbed_info@jbed$readBed(
					.jarray(snp_vec))
			mat_ref = geno_data_ref$getGeno()
#			class(mat_ref) = c("jarrayRef", "jobjRef")
			res = getJArray(mat_ref = mat_ref)
			snp_names = geno_data_ref$getSnpNames()
			if(snp_names_as_colnames) {
				res = setNames(res, snp_names)
			}
			if(fid_iid) {
				res = cbind(fidIid(rbed_info@pl_info), res)
				res$FID = as.character(res$FID)
				res$IID = as.character(res$IID)
			}
			res
		})





#' @rdname readBed
#' @export
setMethod("readBed",
		signature(rbed_info = "RbedInfoC", snp_vec = "missing",
				fid_iid = "missing", snp_names_as_colnames = "missing"),
		function(rbed_info, snp_vec,
				fid_iid, snp_names_as_colnames
		) {
			nsnp = rbed_info@jbed$getnSNPs()
			readBed(rbed_info, 1:nsnp, TRUE, TRUE)
		})

#' @rdname readBed
#' @export
setMethod("readBed",
		signature(rbed_info = "RbedInfoC", snp_vec = "ANY",
				fid_iid = "missing", snp_names_as_colnames = "missing"),
		function(rbed_info, snp_vec,
				fid_iid, snp_names_as_colnames
		) {
			readBed(rbed_info, snp_vec, TRUE, TRUE)
		})

#' @rdname readBed
#' @export
setMethod("readBed",
		signature(rbed_info = "RbedInfoC", snp_vec = "missing",
				fid_iid = "logical", snp_names_as_colnames = "missing"),
		function(rbed_info, snp_vec,
				fid_iid, snp_names_as_colnames
		) {
			nsnp = rbed_info@jbed$getnSNPs()
			readBed(rbed_info,
					snp_vec = 1:nsnp,
					fid_iid = fid_iid,
					snp_names_as_colnames = TRUE)
		})

#' @rdname readBed
#' @export
setMethod("readBed",
		signature(rbed_info = "RbedInfoC", snp_vec = "ANY",
				fid_iid = "logical", snp_names_as_colnames = "missing"),
		function(rbed_info, snp_vec,
				fid_iid, snp_names_as_colnames
		) {
			readBed(rbed_info, snp_vec, fid_iid, TRUE)
		})

#' @rdname readBed
#' @export
setMethod("readBed",
		signature(rbed_info = "RbedInfoC", snp_vec = "missing",
				fid_iid = "missing", snp_names_as_colnames = "logical"),
		function(rbed_info, snp_vec,
				fid_iid, snp_names_as_colnames
		) {
			nsnp = rbed_info@jbed$getnSNPs()
			readBed(rbed_info, 1:nsnp, TRUE, snp_names_as_colnames)
		})


#' @rdname readBed
#' @export
setMethod("readBed",
		signature(rbed_info = "RbedInfoC", snp_vec = "ANY",
				fid_iid = "missing", snp_names_as_colnames = "logical"),
		function(rbed_info, snp_vec,
				fid_iid, snp_names_as_colnames
		) {
			readBed(rbed_info, snp_vec, TRUE, snp_names_as_colnames)
		})

#' @rdname readBed
#' @export
setMethod("readBed",
		signature(rbed_info = "RbedInfoC", snp_vec = "missing",
				fid_iid = "logical", snp_names_as_colnames = "logical"),
		function(rbed_info, snp_vec,
				fid_iid, snp_names_as_colnames
		) {
			nsnp = rbed_info@jbed$getnSNPs()
			readBed(rbed_info, 1:nsnp, fid_iid, snp_names_as_colnames)
		})

#' Add a "shift" suffix to a stem
#'
#' @param stem character.
#' @param n_shift numeric.
#' @return character.
#' @examples
#' \dontrun{
#' # add suffix to stem
#' shiftedStem("a", 100) == "a_shift_0100"
#' shiftedStem("home/a", 100) == "home/a_shift_0100"
#' shiftedStem("/home/a", 100) == "/home/a_shift_0100"
#' shiftedStem(c("/home/a", "/home/b"), 100) == c("/home/a_shift_0100",
#' 		"/home/b_shift_0100")
#' }
#'
#' @author Kaiyin Zhong, Fan Liu
#' @export
shiftedStem = function(stem, n_shift) {
	sprintf("%s_shift_%04d", stem, n_shift)
}



#' Shift bed files
#'
#' Generates collapsed genotypes by shifting the bed file
#' (i.e. SNP1 collapsed with SNP2, SNP2 collapsed with SNP3, etc,
#' when \code{n_shift == 1}).
#'
#' @param rbed_info RbedInfoC object
#' @param n_shift integer.
#' @param collapse_matrix matrix of integers. See details.
#' @param db_setup logical. Whether to setup SQLite database for .bim, .fam and .frq files.
#'
#' @details Collapsing matrix.
#' The collapse_matrix parameter allows collapsing of two genotypes in
#' a arbitrary way. Each genotype is represented by either 0, 1, 2, or 3:
#' \describe{
#' \item{0}{Homozygote of the minor allele.}
#' \item{1}{NA}
#' \item{2}{Heterozygote.}
#' \item{3}{Homozygote of the major allele.}
#' }
#'
#' The collapsing function is implemented as a matrix lookup function, i.e.
#' \eqn{Collapse(S1, S2) = CollapseMatrix[S1][S2]}.
#'
#' The default collapsing matrix is:
#'
#' \tabular{rrrr}{
#'   0 \tab 0 \tab 0 \tab 0\cr
#'   0 \tab 1 \tab 1 \tab 1\cr
#'   0 \tab 1 \tab 0 \tab 3\cr
#'   0 \tab 1 \tab 3 \tab 3
#' }
#'
#' @return RbedInfoC object, with the shifted bed file path in it.
#'
#' @author Kaiyin Zhong, Fan Liu
#' @export
shiftBed = function(rbed_info, n_shift, db_setup = FALSE, collapse_matrix = NULL) {
	stopifnot(isS4Class(rbed_info, "RbedInfoC"))
	n_shift = as.integer(n_shift)
	if(!is.null(collapse_matrix)) {
		if(! "jrectRef" %in% class(collapse_matrix)) {
			stopifnot(is.integer(collapse_matrix) &&
							all(dim(collapse_matrix) == c(4, 4)))
			collapse_matrix = .jarray(collapse_matrix, dispatch = TRUE)
		}
		rbed_info@jbed$shift(n_shift, collapse_matrix)
	} else {
		rbed_info@jbed$shift(n_shift)
	}
	invisible(
			rbedInfo(bedstem = shiftedStem(rbed_info@pl_info@plink_stem, n_shift),
					db_setup = db_setup)
	)
}

#' Remove files by matching the starting part
#'
#' If \code{x} is a string, then this function matches \code{x*} by globbing.
#' If \code{x} is a "PlInfoC" object, it matches \code{x@@plink_stem*},
#' If \code{x} is a "RbedInfoC" object, it matches \code{x@@pl_info@@plink_stem*}.
#' Otherwise nothing is removed.
#'
#'
#' @param x character, PlInfoC, or RbedInfoC object.
#'
#' @author Kaiyin Zhong, Fan Liu
#' @export
rmFilesByStem = function(x) {
	if(is.character(x)) {
		file.remove(Sys.glob(paste(x, "*", sep = "")))
	} else if(isS4Class(x, "PlInfoC")) {
		file.remove(Sys.glob(paste(x@plink_stem, "*", sep = "")))
	} else if(isS4Class(x, "RbedInfoC")) {
		file.remove(Sys.glob(paste(x@pl_info@plink_stem, "*", sep = "")))
	} else {
		NULL
	}
}

#' Create GCDH task directories by tag
#'
#' The task folder is a subfolder of the value of \code{collenv$.collapsabel_gcdh}.
#' It will be created if it does not yet exist.
#'
#' @param gcdh_tag character. Tag for GCDH task.
#' @return character. Directory of the task.
#'
#' @author Kaiyin Zhong, Fan Liu
#' @export
gcdhDir = function(gcdh_tag) {
	stopifnot(is.character(gcdh_tag) && length(gcdh_tag) == 1)
	d = file.path(collenv$.collapsabel_gcdh, gcdh_tag)
	dir.create2(d)
	d
}






minimalPIndices = function(p) {
	stopifnot(is.matrix(p))
	p[is.na(p)] = 99
	max.col(-1 * p, "first")
}

secondSnpIndices = function(p) {
	minimalPIndices(p) + 0:(nrow(p) - 1)
}

secondSnpIndices1 = function(x) {
	x + 0:(length(x) -1)
}


#' Generate a report from a GCDH run
#'
#' For each p-value from a GCDH run, search for indices of the corresponding SNP pair.
#' Combine statistics from single-SNP approach with GCDH statistics.
#'
#' @param run_res Result from \code{runGcdh}
#' @return path to SQLite database
#'
#' @author Kaiyin Zhong, Fan Liu
#' @export
gcdhReport = function(run_res) {
	stopifnot(all(collenv$.linear_header_default %in% names(run_res)))
	# minimal p and number of tests
	gcdh_p = biganalytics::apply(run_res$P, 1, function(i) {
				min(na.omit(i))
			})
	gcdh_ntests = biganalytics::apply(run_res$P, 1, function(i) {length(na.omit(i))})
	pl_gwas = setupRbed(run_res$pl_gwas)
	maf = getQuery(sqliteFilePl(pl_gwas@pl_info), "select maf from frq order by rowid")
	# chr1, bp1, snp1, maf1, nmiss1, beta1, stat1, p1
	basic_info1 = setNames(cbind(run_res$chr_bp, maf), c("SNP1", "CHR1", "BP1", "MAF1"))
	stats1 = data.frame(
			NMISS1 = run_res$NMISS[, 1],
			BETA1 = run_res$BETA[, 1],
			STAT1 = run_res$STAT[, 1],
			P1 = run_res$P[, 1]
	)
	# chr2, bp2, snp2, maf2, nmiss2, beta2, stat2, p2
	minimal_p_indices = minimalPIndices(run_res$P[,])
	second_snp_indices = secondSnpIndices1(minimal_p_indices)
	basic_info2 = setNames(basic_info1[second_snp_indices, ], c("SNP2", "CHR2", "BP2", "MAF2"))
	stats2 = data.frame(
			NMISS2 = stats1$NMISS1[second_snp_indices],
			BETA2  = stats1$BETA1[second_snp_indices],
			STAT2  = stats1$STAT1[second_snp_indices],
			P2     = stats1$P1[second_snp_indices]
	)
	# gcdh_nmiss, gcdh_beta, gcdh_stat, gcdh_p
	idx = cbind(1:nrow(run_res$NMISS), minimal_p_indices)
	gcdh_nmiss = run_res$NMISS[,][idx]
	gcdh_beta = run_res$BETA[,][idx]
	gcdh_stat = run_res$STAT[,][idx]
	gcdh_p1 = run_res$P[,][idx]
	stopifnot(all(gcdh_p1 == gcdh_p))
	# combined all stats in one data.frame
	res = cbind(basic_info1, stats1, basic_info2, stats2,
			data.frame(NMISS = gcdh_nmiss,
					BETA = gcdh_beta,
					STAT = gcdh_stat,
					P = gcdh_p,
					NTEST = gcdh_ntests))
	# write report to SQLite database
	gcdh_report_sqlite_db = sqliteFileGcdh(pl_gwas@gwas_tag, "gcdh_report")
	db = RSQLite::dbConnect(
			RSQLite::SQLite(),
			gcdh_report_sqlite_db
	)
	tryCatch({
				RSQLite::dbWriteTable(db, "gcdh_report", res, overwrite = TRUE)
			}, finally = {
				RSQLite::dbDisconnect(db)
			})
	rm(res)
	gcdh_report_sqlite_db
}





#' Filter a PlGwasC object by the results of a \code{plink --assoc} run
#'
#' This is meant for reduction in computational burden. The \code{plink --assoc} does not
#' accept covariates makes some assumptions accordingly, and thus runs faster than \code{--linear} and
#' \code{--logistic}. SNPs that does not produce a p-value more significant than a user-set threshold will
#' be filtered out. A new PLINK file is made and a corresponding new PlGwasC object is returned.
#'
#' @param pl_gwas PlGwasC object
#' @param plink_out_stem character. Output plink file stem (without .bed extension). The default is to add a "_filtered_{RANDOM_ID}" suffix to the original.
#' @param p_threshold numeric. P-value threshold.
#' @param db_setup logical. Whether to setup the PlGwasC object.
#' @param force logical. Overwrite existing PLINK files.
#' @return a new PlGwasC object. 
#' 
#' @examples 
#' \dontrun{
#' rbed_info = rbedInfo(bedstem = "mmp13", db_setup = FALSE)
#' pl_gwas = plGwas(rbed_info, 
#' 		pheno = "mmp13.phe",
#' 		pheno_name = "Page", 
#' 		gwas_tag = "mmp13_page_sex_age")
#' runGwas(pl_gwas)
#' x = readGwasOut(pl_gwas, c("SNP", "P"), rmGwasOut = FALSE)
#' pl_gwas1 = assocFilter(pl_gwas, p_threshold = 0.001)
#' runGwas(pl_gwas1)
#' x1 = readGwasOut(pl_gwas1, c("SNP", "P"), rmGwasOut = FALSE)
#' y = dplyr::inner_join(x, x1, by = "SNP")
#' all(y$P.x == y$P.y)
#' all(y$P.y < 0.001)
#' }
#'
#' @author Kaiyin Zhong, Fan Liu
#' @export
assocFilter = function(pl_gwas, plink_out_stem = NULL, p_threshold = 0.1, db_setup = FALSE, force = TRUE) {
	stopifnot(is.numeric(p_threshold) && p_threshold > 0 && p_threshold < 1)
	rand_id = randomString(6)
	if(is.null(plink_out_stem))
		plink_out_stem = sprintf("%s_filtered_%s", pl_gwas@pl_info@plink_stem, rand_id)
	target_bed = sprintf("%s.bed", plink_out_stem)
	if(file.exists(target_bed) && !force) stopFormat("File already exists: ", target_bed)
	pl_gwas1 = chGwasTag(pl_gwas, paste0(pl_gwas@gwas_tag, "_", rand_id))
#	pl_gwas1 = setOptModel(pl_gwas1, "assoc")
	runGwas(pl_gwas1)
	assoc_out = readGwasOut(pl_gwas1, c("SNP", "P"))
	assoc_out = assoc_out[which(assoc_out$P < p_threshold), ]
	if(nrow(assoc_out) == 0) {
		stop("No SNPs left after filtering. p_threshold too stringent?")
	}
	snp_list_file = file.path(gwasDir(pl_gwas1), "assoc_snp_list.txt")
	write.table(assoc_out$SNP, quote = FALSE, col.names = FALSE, row.names = FALSE, file = snp_list_file)
	plinkr(bfile = pl_gwas1@pl_info@plink_stem, extract = snp_list_file, make_bed = "", out = plink_out_stem)
	rbed_info = rbedInfo(plink_out_stem, db_setup)
	removeTag(pl_gwas1@gwas_tag)
	pl_gwas2 = plGwas(rbed_info,
			pheno = pl_gwas@opts$pheno,
			pheno_name = pl_gwas@opts$pheno_name,
			covar_name = pl_gwas@opts$covar_name %||% "",
			gwas_tag = pl_gwas@gwas_tag)
	pl_gwas2
}


#' Run GCDH analysis
#'
#' Runs GCDH over the given PlGwasC object.
#' The PlGwasC object is first filtered by p-values from a \code{plink --assoc} run if a p-value threshold is given.
#' New PlGwasC objects are generated by shifting the PLINK bed file (e.g. shift1.bed, shift2.bed, ...) one by one.
#' A GWAS is run for each of these PlGwasC objects and results are collected into big.matrix files.
#'
#' @param pl_gwas PlGwasC object
#' @param n_shift integer. Maximum shift number.
#' @param gwas_col_select character. Columns to read from a GWAS output file. Default to \code{collenv$.linear_header_default}
#' @param collapse_matrix matrix. 4 by 4 matrix used for generating collapsed genotypes.
#' @param rm_shifted_files logical. Whether to remove shifted bed files after analysis is done.
#' @param dist_threshold integer. SNPs beyond this distance will be ignored. Default to 500kb.
#' @return A list with the following members: (1) the input PlGwasC object. (2) an info data frame with CHR, BP and SNP columns. (3) One big.matrix object for each of the names in \code{gwas_col_select}
#' @importFrom bigmemory big.matrix attach.big.matrix filebacked.big.matrix
#' @importFrom biganalytics apply colmin
#' @author Kaiyin Zhong, Fan Liu
#' @export
runGcdh = function(
		pl_gwas,
		n_shift,
		gwas_col_select = NULL,
		collapse_matrix = NULL,
		rm_shifted_files = TRUE,
		dist_threshold = 500e3) {
	if(is.null(gwas_col_select)) {
		gwas_col_select = collenv$.linear_header_default
	}
	# a random suffix makes multiple R session conflicts impossible
	gcdh_tag = sprintf("%s_%s", pl_gwas@gwas_tag, randomString(6))
	pl_gwas = chGwasTag(pl_gwas, gcdh_tag)
	# run the initial GWAS and store results in first column
	runGwas(pl_gwas)
	gwas_out = readGwasOut(pl_gwas, gwas_col_select)
	pl_gwas = setupRbed(pl_gwas)
	nsnps = nSnpPl(pl_gwas@pl_info)
	for(col_name in gwas_col_select) {
		assign(
				col_name,
				gcdhBmCreate(gcdh_tag, col_name, nsnps)
		)
		do.call("[<-", list(x = get(col_name), j = 1, value = gwas_out[, col_name]))
		# reload big matrices after modification
		assign(col_name, bigmemory::attach.big.matrix(bmFilepath(gcdh_tag, col_name, "desc")))
	}
	# run GWAS on shifted bed files
	for(shift_i in 1:n_shift) {
		rbed_info_shifted = shiftBed(pl_gwas, shift_i, FALSE, collapse_matrix)
		shifted_tag = paste(gcdh_tag, "_shift_", shift_i, sep = "")
		pl_gwas_shifted = plGwas(rbed_info_shifted,
				pheno = pl_gwas@opts$pheno,
				pheno_name = pl_gwas@opts$pheno_name,
				covar_name = {
					if("covar_name" %in% names(pl_gwas@opts))
						pl_gwas@opts$covar_name
					else
						""
				},
				gwas_tag = shifted_tag
		)
		runGwas(pl_gwas_shifted)
		gwas_out_shifted = readGwasOut(pl_gwas_shifted, gwas_col_select)
		for(col_name in gwas_col_select) {
			bmAddCol(bmFilepath(gcdh_tag, col_name, "bin"), gwas_out_shifted[, col_name])
			# reload big matrices after modification
			assign(col_name, bigmemory::attach.big.matrix(bmFilepath(gcdh_tag, col_name, "desc")))
		}
		# remove shifted files when requested
		if(rm_shifted_files) {
			file.remove(
					Sys.glob(
							paste(
									rbed_info_shifted@pl_info@plink_stem,
									"*",
									sep = ""
							)
					)
			)
		}
		removeTag(shifted_tag, "gwas")
	}
	# reload big matrices after modification
	for(col_name in gwas_col_select) {
		assign(col_name, bigmemory::attach.big.matrix(bmFilepath(gcdh_tag, col_name, "desc")))
	}
	# mark SNPs that are not in the same chr as NA
	chr_bp = getQuery(sqliteFilePl(pl_gwas@pl_info), "select snp, chr, bp from bim order by rowid")
	for(shift_i in 1:n_shift) {
		dist_i = lagDistance(chr_bp[, "BP"] , shift_i)
		idx = dist_i > dist_threshold
		idx[is.na(idx)] = TRUE
		idx = which(idx)
		chr_diff_i = lagDistance(chr_bp[, "CHR"], shift_i)
		idx1 = chr_diff_i != 0
		idx1[is.na(idx1)] = TRUE
		idx1 = which(idx1)
		idx = unique(c(idx, idx1))
		if(length(idx) > 0) {
			for(col_name in gwas_col_select) {
				do.call("[<-", list(get(col_name), i = idx, j = shift_i + 1, value = NA))
			}
		}
	}
	# make sure there is at least one non-NA p value in each row
	if("P" %in% gwas_col_select) {
		P1 = P[, 1]
		na_idx_P1 = is.na(P1)
		P[na_idx_P1, 1] = 1
	}
	# return stats
	res = list(
			pl_gwas = pl_gwas,
			chr_bp = chr_bp,
			tag = gcdh_tag
	)
	for(col_name in gwas_col_select) {
		res = do.call("[[<-", list(res, col_name, get(col_name)))
	}
	# remove the 0 shift GWAS
	removeTag(pl_gwas@gwas_tag, "gwas")
	res
}



#' Run GCDH over a region
#'
#' A region around some SNP is extracted and GCDH analysis is conducted over that region.
#'
#' @param pl_gwas PlGwasC object
#' @param n_shift integer. Maximum shift number.
#' @param snp character. SNP name
#' @param window numeric. All variants with physical position no more than half the specified kb distance (decimal permitted) from the named variant are loaded.
#' @param out character. Path to the regional bed file (without .bed extension).
#' @param gwas_col_select character. See \code{runGcdh}
#' @param collapse_matrix See \code{runGcdh}
#' @param rm_shifted_files See \code{runGcdh}
#' @param dist_threshold  See \code{runGcdh}
#' @return  See \code{runGcdh}
#'
#' @author Kaiyin Zhong, Fan Liu
#' @export
gcdhRegion = function(pl_gwas,
		n_shift = NULL,
		snp,
		window = 500,
		out = NULL,
		gwas_col_select = collenv$.linear_header_default,
		collapse_matrix = NULL,
		rm_shifted_files = TRUE,
		dist_threshold = 500e3
) {
	suffix = paste0("_", randomString(6), "_window_", window)
	if(is.null(out)) {
		out = paste0(pl_gwas@pl_info@plink_stem, suffix)
	}
	new_tag = paste0(pl_gwas@gwas_tag, suffix)
	plinkr(bfile = pl_gwas@pl_info@plink_stem,
			keep = pl_gwas@opts$pheno,
			snp = snp,
			window = window,
			make_bed = "",
			out = out)
	rbed_info = rbedInfo(out, TRUE)
	pl_gwas = plGwas(rbed_info,
			pheno = pl_gwas@opts$pheno,
			pheno_name = pl_gwas@opts$pheno_name,
			covar_name = pl_gwas@opts$covar_name %||% "",
			gwas_tag = new_tag
	)
	nsnps = nSnpPl(pl_gwas@pl_info)
	if(n_shift >= nsnps || is.null(n_shift)) {
		n_shift = nsnps - 1
	}
	runGcdh(pl_gwas, n_shift, gwas_col_select, collapse_matrix, rm_shifted_files, dist_threshold)
}


#' Get row number of SNPs from their names
#'
#' @param pl_info PlInfoC object.
#' @param snp_names character. Vector of SNP names.
#' @return integer. Vector of row numbers.
#'
#' @author Kaiyin Zhong, Fan Liu
#' @export
snpRowId = function(pl_info, snp_names) {
	if(!isSetup(pl_info)) {
#        print(str(pl_info))
		setup(pl_info)
	}
	getQuery(sqliteFilePl(pl_info),
			sprintf("select rowid from bim where snp in %s order by rowid", strVectorSQLRepr(snp_names)))[, 1]
}
