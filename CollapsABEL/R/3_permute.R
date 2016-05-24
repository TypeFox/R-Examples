#' Validate a phenotype file
#' 
#' @param phe_file character. Phenotype file.
#' @param ... Passed to read.table
#' @return FALSE when the file is invalid, or a data.frame when it is. 
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
validPhe = function(phe_file, ...) {
	phe = read.table(phe_file, stringsAsFactors = FALSE, header = TRUE, ...)
	if(!all(colnames(phe)[1:2] == c("FID", "IID"))) {
		return(FALSE)
	}
	for(i in 3:ncol(phe)) {
		if(!is.numeric(phe[, i])) {
			return(FALSE)
		}
	}
	phe
}

#' Permute a phenotype file
#' 
#' All columns except FID and IID are permuted.
#' 
#' @param phe_file character. Phenotype file.
#' @param out_file character. Path to permuted phenotype file.
#' @param force logical. When set to TRUE, existing file is overwritten.
#' @param valid logical. Whether to validate the phenotype file first. 
#' @param ... Passed to read.table
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
permutePhe = function(phe_file, out_file, force = FALSE, valid = TRUE, ...) {
	if(file.exists(out_file) && !force) {
		stopFormat("File %s already exists!", out_file)
	}
	if(valid) phe = validPhe(phe_file, ...)
	else phe = read.table(phe_file, stringsAsFactors = FALSE, header = TRUE, ...)
	fidiid = phe[, c("FID", "IID")]
	others = phe[, !(colnames(phe) %in% c("FID", "IID"))]
	others = others[sample(nrow(others)), ]
	phe = cbind(fidiid, others)
	write.table(phe, file = out_file, quote = FALSE, row.names = FALSE)
	invisible(phe)
}

#' Run simulations to control type-I error
#' 
#' Simulate a new phenotype N times and run GCDH with each. 
#' The \code{phe_fun} function is used to generate new phenotype file. 
#' When this function is not given, the phenotype file from the PlGwasC object will be
#' permuted and used as the new phenotype file (permutation analysis). Thus when no \code{phe_fun} 
#' is supplied, this function can be used to survey p-values under the null distribution. 
#' A threshold for Genome-wide significance can be calculated from these p-values by 5% (or
#' any other alpha-level) quantile.
#' 
#' @param pl_gwas PlGwasC object
#' @param n_shift integer. \code{n_shift} for each GCDH run.
#' @param n_simu integer. Number of simulations to run.
#' @param phe_fun function. Used to generate new phenotype file.
#' @param dist_threshold See runGcdh.
#' @param p_threshold numeric or NULL. When it's not NULL, the PlGwasC object is filtered by \code{assocFilter} first.
#' @param collapse_matrix See runGcdh.
#' @param rm_shifted_files  See runGcdh.
#' @return A list with the following members: (1) tag of this simulation, can be used to remove related files. (2) 
#' a list of SNP pairs. If "snp_pair" is a member of the result from \code{phe_fun}, then this list will be non-empty, 
#' otherwise it will be empty. (3) a list of reports from all the GCDH analysis. (4) global minimal p-values of the single-SNP
#' approach. (4) global minimal p-values of GCDH.
#' 
#' @author Kaiyin Zhong, Fan Liu
#' @export
runTypeI = function(
		pl_gwas,
		n_shift, 
		n_simu,
		phe_fun = NULL,
		dist_threshold = 500e3,
		p_threshold = NULL,
		collapse_matrix = NULL,
		rm_shifted_files = TRUE
) {
	stopifnot(is.data.frame(validPhe(pl_gwas@opts$pheno)))
	# to avoid collision between multiple R sessions, add a suffix of random chars
	rand_id = randomString(6)
	type_one_tag = sprintf("%s_%s", pl_gwas@gwas_tag, rand_id)
	pl_gwas = chGwasTag(pl_gwas, type_one_tag)
	pl_gwas = setupRbed(pl_gwas)
	nsnps = nSnpPl(pl_gwas@pl_info)
	snp_pairs = list()
	gcdh_reports = list()
	typeI_dir = tag2Dir(type_one_tag, "typeI")
	dir.create2(typeI_dir)
	new_pheno = file.path(typeI_dir, "permuted.phe")
	if(is.null(phe_fun)) {
		phe_fun = function(filename) {
			permutePhe(pl_gwas@opts$pheno, filename, TRUE, FALSE)
		}
	}
	for(simu_i in 1:n_simu) {
		phe_fun_ret = phe_fun(new_pheno) 
		if("snp_pair" %in% names(phe_fun_ret)) {
			snp_pairs[[length(snp_pairs) + 1]] = phe_fun_ret$snp_pair
		}
		new_pl_gwas = pl_gwas
		new_pl_gwas@opts$pheno = new_pheno
		new_pl_gwas = chGwasTag(new_pl_gwas, paste0(new_pl_gwas@gwas_tag, "_simulation_", simu_i))
		if(!is.null(p_threshold)) {
			new_pl_gwas = assocFilter(new_pl_gwas, p_threshold = p_threshold)
		}
		gcdh_res = runGcdh(
				pl_gwas = new_pl_gwas, 
				n_shift = n_shift, 
				collapse_matrix = collapse_matrix, 
				rm_shifted_files = rm_shifted_files, 
				dist_threshold = dist_threshold
		)
		gcdh_reports[[length(gcdh_reports) + 1]] = getQuery(gcdhReport(gcdh_res), "select snp1, snp2, p1, p2, p from gcdh_report order by rowid")			
		removeTag(gcdh_res$tag, "gcdh")
		if(!is.null(p_threshold)) {
			rmPlinkFiles(new_pl_gwas@pl_info@plink_stem)
		}
	}
	pminmins = sapply(gcdh_reports, function(report_i) {
				min(report_i$P, na.rm = T)
			})
	pbasemins = sapply(gcdh_reports, function(report_i) {
				m1 = min(report_i$P1, na.rm = T)
				m2 = min(report_i$P2, na.rm = T)
				min(m1, m2, na.rm = T)
			})
	list(
			tag = type_one_tag,
			snp_pairs = snp_pairs,
			gcdh_reports = gcdh_reports,
			pbasemins = pbasemins,
			pminmins = pminmins
	)
}


dummyPhe = function(rbed_info, filename) {
	phe = fidIid(rbed_info@pl_info)
	phe$y = rnorm(nrow(phe))
	write.table(phe, quote = FALSE, row.names = FALSE, file = filename)	
}

dummyTypeI <- function(rbed_info, 
		n_shift, 
		n_simu, 
		...) {
	phe_fun = function(filename) dummyPhe(rbed_info, filename)
	dummy_tag = paste0("dummy_", randomString(6))
	dummy_dir = tag2Dir(dummy_tag, "dummy")
	dir.create2(dummy_dir)
	phe_file = file.path(dummy_dir, paste0("phe_", randomString(6), ".txt"))
	phe_fun(phe_file)
	pl_gwas = plGwas(pl_gwas = rbed_info, 
			pheno = phe_file,
			pheno_name = "y", 
			gwas_tag = dummy_tag)
	typeI_res = runTypeI(pl_gwas = pl_gwas, 
			n_shift = n_shift, 
			n_simu = n_simu, 
			phe_fun = phe_fun, 
			...)
	unlink(dummy_dir, recursive = TRUE)
	typeI_res
}


#' GCDH power analysis
#' 
#' This function makes use of \code{runTypeI}. 
#' Random phenotypes are used to survey p-values under the null hypothesis (SNPs are not associated phenotype), 
#' and genome-wide significance thresholds for single-SNP approach and GCDH are calculated by a user given 
#' alpha-level.
#' A custom \code{phe_fun} is supplied 
#' for simulating a phenotype associated with a certain pair of SNPs. 
#' Total number of such simulations is set by the n_simu parameter. In each simulation 4 p-values are generated: 
#' 
#' P_single: p-values from single-SNP approach.
#' 
#' P_GCDH:  p-values from GCDH.
#' 
#' P_(single,no causal): p-values from single-SNP approach when causal SNPs are untyped.
#' 
#' P_(GCDH,no causal): p-values from GCDH when causal SNPs are untyped.
#' 
#' When all simulations are finished, 
#' 4 vectors of p-values are obtained: P_single_vec, P_GCDH_vec, P_(single,no causal)_vec, P_(GCDH,no causal)_vec.
#' The power for each of the category (single-SNP, single-SNP without causal genotypes, GCDH, GCDH without causal genotypes)
#' are proportions of these vectors that are more significant than the genome-wide significance thresholds 
#' we have obtained.
#' 
#' @param rbed_info RbedInfoC object
#' @param n_shift integer. \code{n_shift} for each GCDH run.
#' @param n_simu integer. Number of simulations to run.
#' @param maf_min numeric. Lower limit of MAF interval. 
#' @param maf_max numeric. Upper limit of MAF interval. 
#' @param r_limit numeric. Upper limit of correlation coefficient between the two causal SNPs.
#' @param beta numeric. Effect size of the simulated phenotype.
#' @param collapse_matrix See runGcdh.
#' @param dist_threshold See runGcdh.
#' @param alpha_level numeric. Control type-I error rate at this level.
#' 
#' @author Kaiyin Zhong
#' @export
gcdhPower = function(rbed_info, n_shift, n_simu, maf_min, maf_max,
		r_limit, beta, 
		collapse_matrix = NULL, 
		dist_threshold = 5e5,
		alpha_level = 0.05
) {
	rbed_info = setupRbed(rbed_info)
	snp_frq = getQuery(sqliteFilePl(rbed_info@pl_info), 
			sprintf("select chr,snp,bp,maf from bim join frq using (chr,snp) where maf > %s and maf < %s order by chr,bp ", maf_min, maf_max))
	phe_fun1 = function(filename) {
		chooseSnps = function() {
			# choose chr, and extract data of that chr
			chr = sample(snp_frq$CHR, 1)
			snp_frq1 = snp_frq[snp_frq$CHR == chr, ]
			# sample for snp1
			snp1 = snp_frq1$SNP[sample(1:nrow(snp_frq1), 1)]
			snp1_idx = which(snp_frq1$SNP == snp1)
			# further extraction, we only need the SNPs within the window size (n_shift)
			end_idx = min(c(snp1_idx + n_shift, nrow(snp_frq1)))
			if((end_idx - snp1_idx) < 3) { 
				return(NULL) # make sure the pool is big enough
			}
			snp_frq1 = snp_frq1[snp1_idx:end_idx, ]
			# further extraction, we only need the SNPs that are within a certain distance (dist_threshold)
			bp1 = snp_frq1[snp_frq1$SNP == snp1, ]$BP[1]
			bpdiff = snp_frq1$BP - bp1
			bpdiff_idx = which(bpdiff < dist_threshold)
			if(length(bpdiff_idx) < 3) {
				return(NULL)
			}
			snp_frq1 = snp_frq1[bpdiff_idx, ]
			# read genotypes
			## choose 10 SNPs as the pool
			snps = snp_frq1$SNP
			if(length(snps) > 10) {
				snps = sample(snps, 10)
			}
			geno = readBed(rbed_info, 
					snp_vec = snps)
			r_vec = apply(geno[, 3:ncol(geno)], 2, function(g) {
						idx = !is.na(g) & !is.na(geno[, 3])
						cor(g[idx], geno[idx, 3])
					})
			r_idx = which(abs(r_vec) < r_limit) + 2
			if(length(r_idx) == 0) {
				return(NULL)
			}
			snps_after_r_filter = colnames(geno)[r_idx]
			snp_pair = c(snps[1], sample(snps_after_r_filter, 1))
			geno[, c("FID", "IID", snp_pair)]
		}
		chooseSnpsLoop = function() {
			geno_pair = chooseSnps()
			n_try = 1 
			while(is.null(geno_pair) && n_try < 21) {
				geno_pair = chooseSnps()
				message(sprintf("No suitable SNPs found in this round of sampling, retry %d", n_try))
				n_try = n_try + 1
			}
			if(is.null(geno_pair)) {
				stop("No suitable SNPs found. Consider using less stringent parameters.")
			}
			geno_pair
		}
		geno = chooseSnpsLoop()
		geno1 = geno[, 3]
		geno2 = geno[, 4]
		geno$y = ((geno1 + geno2) >= 2) * beta + rnorm(length(geno1))
		write.table(geno[, c(1, 2, 5)], file = filename, 
				quote = FALSE, row.names = FALSE)
		list(snp_pair = colnames(geno)[3:4])
	}
	phe_dummy = tempfile()
	dummyPhe(rbed_info, phe_dummy)
	pl_gwas = plGwas(rbed_info, 
			pheno = phe_dummy, 
			pheno_name = "y", 
			gwas_tag = "power_est")
	typeI_res = runTypeI(pl_gwas, 
			n_shift = n_shift, 
			n_simu = n_simu, 
			phe_fun = phe_fun1, 
			dist_threshold = dist_threshold, 
			collapse_matrix = collapse_matrix
	)
	pmins_no_causal = t(
			sapply(1:length(typeI_res$snp_pairs), function(i) {
						snp_pair = typeI_res$snp_pairs[[i]]
						gcdh_report = typeI_res$gcdh_reports[[i]]
						idx = (! (gcdh_report$SNP1 %in% snp_pair)) & (! (gcdh_report$SNP2 %in% snp_pair))
						gcdh_report = gcdh_report[idx, ]
						res = apply(gcdh_report[, c("P1", "P2", "P")], 2, function(col_iter) {
									min(col_iter, na.rm = TRUE)
								})
						res = c(min(res[c("P1", "P2")]), res["P"])
						names(res) = c("PSINGLE", "P")
						res
					})
	)
	typeI_res$pbasemins_no_causal = pmins_no_causal[,"PSINGLE"]
	typeI_res$pminmins_no_causal = pmins_no_causal[, "P"]
	dummy_res = dummyTypeI(rbed_info, n_shift, n_simu, 
			dist_threshold = dist_threshold, 
			collapse_matrix = collapse_matrix)[c("pbasemins", "pminmins")]
	alpha_single = quantile(dummy_res$pbasemins, alpha_level)
	alpha_gcdh = quantile(dummy_res$pminmins)
	pow_single = mean(typeI_res$pbasemins < alpha_single, na.rm = TRUE)
	pow_gcdh = mean(typeI_res$pminmins < alpha_gcdh, na.rm = TRUE)
	pow_single_nc = mean(typeI_res$pbasemins_no_causal < alpha_single, na.rm = TRUE)
	pow_gcdh_nc = mean(typeI_res$pminmins_no_causal < alpha_gcdh, na.rm = TRUE)
	res = data.frame(pow_single = pow_single, 
			pow_gcdh = pow_gcdh, 
			pow_single_nc = pow_single_nc, 
			pow_gcdh_nc = pow_gcdh_nc)
	res
}

#' Generate phenotype file from a fam file 
#' @param famfile Character. Path of fam file.
#' @param n_components Integer. Number of principle components to generate. 
#' @return Phenotype data.frame. The data frame contains the FID, IID, SEX, AFFECTEDNESS columns of the fam file, plus principle components of genetic information.
#' 
#' @author kaiyin
#' @export
makePhe = function(famfile, n_components) {
	base_filename = tools::file_path_sans_ext(famfile)
	phe_file = paste0(base_filename, ".phe")
	pca_file = paste0(base_filename, ".eigenvec")
	phe = readFam(famfile, c("FID", "IID", "SEX", "PHE"))
	plinkr(bfile = base_filename, pca = n_components, out = base_filename)
	pca = read.table(pca_file, header = FALSE, stringsAsFactors = FALSE)
	colnames(pca) = c("FID", "IID", paste0("pca", 1:n_components))
	phe$FID = as.character(phe$FID)
	phe$IID = as.character(phe$IID)
	pca$FID = as.character(pca$FID)
	pca$IID = as.character(pca$IID)
	phe = dplyr::left_join(phe, pca)
	write.table(phe, file = phe_file, quote = FALSE, row.names = FALSE)
	phe
}

