
.fishersMethod <- function(value, orig_df) stats::pchisq(-2 * value,df=2*orig_df,lower.tail=FALSE)

.normalizeData <- function(orig_data, norm_method) {

  print(paste("normalizing data by", norm_method, sep=" "))
  if (norm_method == "None") {
    subtract = rep(0, nrow(orig_data))
  }
  else if (norm_method == "All") {
    subtract <- colMeans(orig_data[, 2:length(orig_data)])
  }
  else if (norm_method %in% colnames(orig_data)) {
    ind <- which(names(orig_data)==norm_method)
    if (length(ind) > 1)
      subtract <- colMeans(orig_data[, ind])
    else
      subtract <- orig_data[, ind]
  }
  else stop("only All or None or name of individual from the dataset is considered.")

  normed_data <- orig_data - subtract

  return(normed_data)
}

#### QTL analysis using anova ####

.qtlAnalysis <- function(exp_data,
                         genotype_data,
                         genes_to_use=c(),
                         transpose.file=F, calc.max=F) {

  # reading the SNP file

  snp_info <- genotype_data[, 1:5]
  snp_data <- genotype_data[, 6:length(genotype_data)]

  if (length(genes_to_use) > 0 ) {
    exp_data <- exp_data[rownames(exp_data) %in% genes_to_use, ]
  }

  if (transpose.file) {
    exp_data <- as.data.frame(t(exp_data))
  } else {
    exp_data <- as.data.frame(exp_data)
  }

  temp_data_name <- paste("`", colnames(exp_data), "`", sep="")

  print("running QTL analysis")
  log_pValue_res <- matrix(nrow=ncol(exp_data), ncol = nrow(snp_data),
                           dimnames=list(colnames(exp_data), rownames(snp_data)))

  ind_In_SNPs_File <- match(rownames(exp_data), colnames(snp_data),nomatch=0)

  snpData_filtered <- snp_data[, ind_In_SNPs_File]
  snpData_filtered <- as.data.frame(t(snpData_filtered))

  print("going over all possible snps to find the best qtl")
  ## going over all possible snps to find the best qtl
  for (j in c(1:ncol(snpData_filtered))) {
    formula <- stats::as.formula( paste0("cbind(", paste(temp_data_name,collapse=","), ")~snpData_filtered[, j]") )

    pval_sum <- summary(stats::aov(formula, exp_data))
    for (i in c(1:ncol(exp_data)))
      log_pValue_res[i, j] <- -log10(pval_sum[[i]][1,5])
  }

  print("end of qtl analysis")
  res_of_cell <- list()

  if (calc.max) {
    max_of_each_row <- apply(log_pValue_res, 1, max)
    for (i in c(1:nrow(log_pValue_res))) {
      name_of_cell <- names(exp_data)[i]
      high_markers <- row.names(snp_info)[which(log_pValue_res[i,]==max_of_each_row[i])]
      high_markers_fomratted_1 <- paste(high_markers, collapse=", ")
      high_markers_fomratted_2 <- paste("[", high_markers_fomratted_1, "]", sep="")
      chrom <- snp_info[which(log_pValue_res[i,]==max(log_pValue_res[i,]))[1], 2]
      place_in_genome <- snp_info[which(log_pValue_res[i,]==max_of_each_row[i])[1], 3]

      res_of_cell[[name_of_cell]] = c(high_markers_fomratted_2, chrom, place_in_genome, max_of_each_row[i])
    }

  }
  return(list("logPval"=log_pValue_res, "max.results"=res_of_cell))

}

############################

.create_iQTL_association_scores <- function(reference_data,
                                            mix_data,
                                            marker_set,
                                            genotype_data) {

  decon_average <- dcq(reference_data = reference_data, mix_data = mix_data,
                       marker_set = marker_set)$average

  res_mat <- .qtlAnalysis(exp_data = decon_average,
                          genotype_data=genotype_data)$logPval

  res_mat[is.na(res_mat)] <- 0

  return(res_mat)

}

.combine_iQTL_association_scores <- function(combined_results) {
  reg_pval_from_minus_log = function(x) 10^(-x)
  mat_log = lapply(combined_results, reg_pval_from_minus_log)
  mat_ln = lapply(mat_log, log)

  sum_mat <- apply(simplify2array(mat_ln), 1:2, sum)

  fisher_mat <- .fishersMethod(sum_mat, length(combined_results))
  fisher_mat_log <- -log10(fisher_mat)

  return(fisher_mat_log)
}


