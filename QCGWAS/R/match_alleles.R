match_alleles <-
function(dataset, ref_set, HQ_subset,
                          dataname = "dataset", ref_name = "reference",
                          unmatched_data = !all(dataset$MARKER %in% ref_set$SNP),
                          check_strand = FALSE, save_mismatches = TRUE, delete_mismatches = FALSE,
                          delete_diffEAF = FALSE, threshold_diffEAF = 0.15,
                          check_FRQ = TRUE, check_ambiguous = FALSE,
                          plot_FRQ = FALSE, plot_intensity = FALSE,
                          plot_if_threshold = FALSE, threshold_r = 0.95,
                          return_SNPs = FALSE, return_ref_values = FALSE, 
                          header_translations, header_reference,
                          save_name = dataname, save_dir = getwd(), use_log = FALSE, log_SNPall = nrow(dataset)) {
  # Note that log_SNPall is meant exclusively for the log files. It
  # won't necessarily equal the number of SNPs passed to the function.
  
  # Part 1: checking headers & merging data  
  #	note that data_col/ref_col have different contents,
  #	depending on whether check_header/check_ref is TRUE or FALSE
  if(delete_diffEAF & !check_FRQ) stop("cannot remove aberrant allele-frequencies when check_FRQ is FALSE")
  check_header <- !missing(header_translations)
  check_ref <- !missing( header_reference)
  col_effect <- TRUE
  
  
  data_std <- c("MARKER", "EFFECT_ALL", "OTHER_ALL", "STRAND", "EFFECT", "EFF_ALL_FREQ")[c(TRUE, TRUE, TRUE, check_strand, TRUE, check_FRQ)]
  if(check_header) {
    if(any(duplicated(header_translations[ ,2]))) stop("duplicated elements in header_translations, column 2")
    data_h <- toupper(colnames(dataset))
    data_col <- integer(length = length(data_std))
    for(forI in 1:length(data_std)) {
      data_cur <- identify_column(data_std[forI], header_translations, data_h)
      if(length(data_cur) == 1L) {
        data_col[forI] <- data_cur
      } else {
        if(length(data_cur) == 0L) { 
          if(data_std[forI] == "EFFECT") { col_effect <- FALSE
          } else { stop(paste("Cannot identify data column:", data_std[forI])) }
        } else { stop(paste("Multiple data columns identified as:", data_std[forI])) }
      }
    }
    if(!col_effect) {
      data_std <- c("MARKER", "EFFECT_ALL", "OTHER_ALL", "STRAND", "EFFECT", "EFF_ALL_FREQ")[c(TRUE, TRUE, TRUE, check_strand, FALSE, check_FRQ)]
      data_col <- data_col[data_col != 0L]
    }
  } else {
    if(!("EFFECT" %in% colnames(dataset))) {
      col_effect <- FALSE
      data_std <- c("MARKER", "EFFECT_ALL", "OTHER_ALL", "STRAND", "EFFECT", "EFF_ALL_FREQ")[c(TRUE, TRUE, TRUE, check_strand, FALSE, check_FRQ)]
    }
    data_col <- data_std
  }
  
  ref_std <- c("SNP", "MINOR", "MAJOR", "MAF")[c(TRUE, TRUE, TRUE, check_FRQ)]
  if(check_ref) {
    if(any(duplicated(header_reference[ ,2]))) stop("duplicated elements in header_reference, column 2")
    ref_h	 <- toupper(colnames(ref_set))
    ref_col <- integer(length = length(ref_std))
    for(forI in 1:length(ref_std)) {
      ref_cur <- identify_column(ref_std[forI], header_reference, ref_h)
      if(length(ref_cur) == 1L) {
        ref_col[forI] <- ref_cur
      } else {
        if(length(ref_cur) == 0L) { stop(paste("Cannot identify reference column:", ref_std[forI]))
        } else {			    stop(paste("Multiple reference columns identified as:", ref_std[forI])) }
      }
    }
  } else { ref_col <- ref_std }
  
  if(any(is.na(dataset[ , data_col[1]]))) { stop("missing SNP IDs in marker name column") }
  if(any(is.na(ref_set[ ,	ref_col[1]]))) { stop("missing SNP IDs in reference SNP name column") }
  if(any(duplicated(dataset[ , data_col[1]]))) { stop("duplicate SNPs in marker name column") }
  if(any(duplicated(ref_set[ , ref_col[1]]))) { stop("duplicate SNPs in reference SNP name column") }
  
  # The awkward temp_order column is necessary because merge does not maintain the original order of
  # x when there are SNPs that have no match in y.
  order_col <- length(data_std) + 1L
  dataset <- merge(x = cbind(dataset[ , data_col], 1:nrow(dataset)), y = ref_set[ , ref_col], by.x = 1, by.y = 1, all.x = TRUE, all.y = FALSE, sort = FALSE)
  if(check_header | check_ref) { colnames(dataset) <- c(data_std, "temp_order", ref_std[-1])
  } else {				 colnames(dataset)[order_col] <- "temp_order" }
  if(unmatched_data) { dataset <- dataset[order(dataset$temp_order), ] }
  
  missing_list <- is.na(dataset$EFFECT_ALL) | is.na(dataset$OTHER_ALL) | is.na(dataset$MINOR) | is.na(dataset$MAJOR)
  SNPn_missing <- sum(missing_list)
  if(SNPn_missing > 0L) {
    SNPn_missing_data <- sum(is.na(dataset$EFFECT_ALL) | is.na(dataset$OTHER_ALL))
    SNPn_missing_ref	<- sum(is.na(dataset$MINOR	) | is.na(dataset$MAJOR	))
    if(SNPn_missing_data > 0L) {
      if(use_log) { save_log(phaseL = 3L, checkL = paste("allele match -", ref_name), typeL = "alleles", SNPL = SNPn_missing_data, allSNPs = log_SNPall, actionL = "excluded from test", noteL = "no allele information: cannot match these SNPs", fileL = paste(save_dir, save_name, sep = "/")) }
      print(paste(" - - missing allele-data in", SNPn_missing_data, "SNPs"), quote = FALSE)
    }
    if(SNPn_missing_ref > 0L) {
      if(use_log) { save_log(phaseL = 3L, checkL = paste("allele match -", ref_name), typeL = "incomplete reference", SNPL = SNPn_missing_ref, allSNPs = log_SNPall, actionL = "excluded from test", noteL = "reference does not contain alleles for these SNP", fileL = paste(save_dir, save_name, sep = "/")) }
      print(paste(" - - incomplete reference: no alleles found for", SNPn_missing_ref, "SNPs"), quote = FALSE)
    }
  } else {
    SNPn_missing_data <- 0L
    SNPn_missing_ref	<- 0L
  }
  
  # Part 2: checking strand column for negative-strand SNPs
  if(check_strand) {
    min_list <- dataset$STRAND == "-" & !is.na(dataset$STRAND)
    SNPn_min <- sum(min_list)
    SNPn_min_SS <- 0L
    SNPn_min_MM <- 0L
    if(SNPn_min > 0L) { dataset[min_list, c("EFFECT_ALL", "OTHER_ALL", "STRAND")] <- switch_strand(dataset[min_list, c("EFFECT_ALL", "OTHER_ALL", "STRAND")], strand_col = TRUE) }
  } else {
    SNPn_min <- NA
    SNPn_min_SS <- NA
    SNPn_min_MM <- NA
  }
  
  # Part 3: checking alleles for strand-mismatches
  switch_list <- !( ( dataset$EFFECT_ALL == dataset$MINOR & dataset$OTHER_ALL == dataset$MAJOR ) |
                      ( dataset$OTHER_ALL == dataset$MINOR & dataset$EFFECT_ALL == dataset$MAJOR ) | missing_list )
  SNPn_switch <- sum(switch_list)
  if(SNPn_switch > 0L) {
    if(check_strand & SNPn_min > 0L) { SNPn_min_SS <- sum(min_list[switch_list]) }
    dataset[switch_list, c("EFFECT_ALL", "OTHER_ALL")] <-
      switch_strand(dataset[switch_list, c("EFFECT_ALL", "OTHER_ALL")], strand_col = FALSE)
    
    # Part 3b: checking for mismatching SNPs (i.e. SNPs whose alleles do not match even after strand-switch)
    #	Depending on the settings, these are switched back to how they were before entering part 3,
    #	saved to a file, and deleted.
    mismatch_list <- !( ( dataset$EFFECT_ALL == dataset$MINOR & dataset$OTHER_ALL == dataset$MAJOR ) |
                          ( dataset$OTHER_ALL == dataset$MINOR & dataset$EFFECT_ALL == dataset$MAJOR ) | missing_list )
    SNPn_mismatch <- sum(mismatch_list)
    
    if(SNPn_mismatch > 0L) {
      missing_list[mismatch_list] <- TRUE
      if(save_mismatches | !delete_mismatches) {
        dataset[mismatch_list, c("EFFECT_ALL", "OTHER_ALL")] <-
          switch_strand(dataset[mismatch_list, c("EFFECT_ALL", "OTHER_ALL")], strand_col = FALSE)
      }
      if(save_mismatches) {
        if(check_strand) {
          negative_strand_correction <- min_list[mismatch_list]
          if(SNPn_min > 0L) { SNPn_min_MM <- sum(negative_strand_correction) }
        }
        write.table( dataset[mismatch_list, -order_col ], paste0(save_dir, "/", save_name, "_SNPs_mismatches-", ref_name, ".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
      } else if(check_strand & SNPn_min > 0L) { SNPn_min_MM <- sum(min_list[mismatch_list]) }
      if(delete_mismatches) { dataset[mismatch_list, "EFFECT_ALL"] <- NA }
      if(use_log) { save_log(phaseL = 3L, checkL = paste("allele match -", ref_name), typeL = "allele mismatch", SNPL = SNPn_mismatch, allSNPs = log_SNPall, actionL = if(delete_mismatches) "SNPs removed" else "-", noteL = "Allele mismatch with reference; strand correction did not correct this", fileL = paste(save_dir, save_name, sep = "/"))
      } else { print(paste(" - - warning:", SNPn_mismatch, "mismatches found - strand correction did not correct this"), quote = FALSE) }
    }
  } else { SNPn_mismatch <- 0L }
  
  # Part 4: matching allele configuration with reference
  flip_list <- dataset$OTHER_ALL == dataset$MINOR & dataset$EFFECT_ALL == dataset$MAJOR & !missing_list
  SNPn_flip <- sum(flip_list)
  if(SNPn_flip > 0L) {
    temp_al2 <- dataset$OTHER_ALL
    dataset$OTHER_ALL <- ifelse(flip_list, dataset$EFFECT_ALL, dataset$OTHER_ALL)
    dataset$EFFECT_ALL <- ifelse(flip_list, temp_al2, dataset$EFFECT_ALL)
    if(col_effect) { dataset$EFFECT <- ifelse(flip_list, -dataset$EFFECT, dataset$EFFECT) }
    if(check_FRQ) { dataset$EFF_ALL_FREQ <- ifelse(flip_list, 1 - dataset$EFF_ALL_FREQ, dataset$EFF_ALL_FREQ) }
  }
  
  # Part 5: checking for "ambiguous" SNPs, i.e. SNPs with an A/T or C/G allele configuration
  ambiguous_list <- (( dataset$EFFECT_ALL == "A" & dataset$OTHER_ALL == "T" ) |
                       ( dataset$EFFECT_ALL == "T" & dataset$OTHER_ALL == "A" ) |
                       ( dataset$EFFECT_ALL == "C" & dataset$OTHER_ALL == "G" ) |
                       ( dataset$EFFECT_ALL == "G" & dataset$OTHER_ALL == "C" ) ) & !missing_list
  SNPn_ambiguous <- sum(ambiguous_list)
  if(check_FRQ & SNPn_ambiguous > 0L) { SNPn_suspect <- sum(ambiguous_list & !is.na(dataset$EFF_ALL_FREQ) & !is.na(dataset$MAF) & ( ( dataset$EFF_ALL_FREQ > 0.65 & dataset$MAF < 0.35 ) | ( dataset$EFF_ALL_FREQ < 0.35 & dataset$MAF > 0.65 ) ))
  } else { SNPn_suspect <- if(check_FRQ) 0L else NA }
  
  # Part 6: correlating allele frequency with reference
  
  if(check_FRQ) {
    FRQ_all <- dataset[ , "EFF_ALL_FREQ"]
    if(SNPn_missing + SNPn_mismatch > 0L) { FRQ_all[missing_list] <- NA }
    diffEAF_list <- abs(FRQ_all - dataset$MAF) > threshold_diffEAF & !is.na(FRQ_all) & !is.na(dataset$MAF)
    SNPn_diffEAF <- sum(diffEAF_list)
    if(SNPn_diffEAF > 0L) {
      if(use_log) { save_log(phaseL = 3L, checkL = paste("allele match -", ref_name), typeL = "allele frequency", SNPL = SNPn_diffEAF, allSNPs = log_SNPall, actionL = if(delete_diffEAF) "Markers removed" else "-", noteL = paste("Allele frequencies deviate from those in", ref_name,"( > ", 100 * threshold_diffEAF, "% )"), fileL = paste(save_dir, save_name, sep = "/"))
      } else { print(paste(" - - warning:", SNPn_diffEAF, "SNPs whose allele-frequency differs strongly from reference"), quote = FALSE) }
      if(delete_diffEAF) {
        write.table(data.frame(dataset[diffEAF_list, -order_col], DIFFERENCE = abs(dataset$EFF_ALL_FREQ[diffEAF_list] - dataset$MAF[diffEAF_list]) ), paste0(save_dir, "/", save_name, "_SNPs_EAFdifferent-", ref_name, ".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
        dataset[diffEAF_list, "EFFECT_ALL"] <- NA
      }
    }
    
    # Creating/editing HQ list
    if(plot_FRQ) {
      if(missing(HQ_subset)) {
        HQ_subset <- !logical(length = nrow(dataset))
        plot_HQ <- FALSE
      } else {
        if(is.numeric(HQ_subset)) {
          if(length(HQ_subset) > nrow(dataset)) stop("HQ_subset cannot be longer than the dataset")
          temp <- logical(length = nrow(dataset))
          temp[HQ_subset] <- TRUE
          HQ_subset <- temp
        } else {
          if(!is.logical(HQ_subset)) stop("HQ_subset is not a logical or numeric vector")
          if(length(HQ_subset) != nrow(dataset)) stop("HQ_subset is not of equal length to the dataset")
        }
        plot_HQ <- TRUE
      }
    }
    
    # Creating correlation/plot function
    plot_allele_cor <- function(dat_FRQ, ref_FRQ, HQ_subset,
                                plot_FRQ, plot_if_threshold, threshold_r,
                                plot_intensity, plot_HQ,
                                name_dir, name_file, name_ref, name_data,
                                ambiguous = NA, 
                                use_log, relN, allN){
      
      # Creating text messages
      if(is.na(ambiguous)) {
        plot_file <- paste0(name_dir, "/", name_file, "_graph_EAF-", name_ref, ".png")
        note_type <- "allele frequency"
        note_correlation1 <- "Allele frequencies correlate poorly with those in"
        note_correlation2 <- " - - warning: allele frequency correlates poorly with reference (r = "
        note_failure <- "insufficient non-missing allele frequencies to correlate SNPs allele frequencies with reference"
      } else{
        if(ambiguous) {
          plot_file <- paste0(name_dir, "/", name_file, "_graph_EAF-", name_ref, "-ambiguous.png")
          note_type <- "ambiguous SNPs"
          note_correlation1 <- "Ambiguous SNPs allele frequencies correlate poorly with those in"
          note_correlation2 <- " - - warning: allele frequency of ambiguous SNPs correlates poorly with reference (r = "
          note_failure <- "insufficient non-missing allele frequencies to correlate ambiguous SNPs allele frequencies with reference"
        } else {
          plot_file <- paste0(name_dir, "/", name_file, "_graph_EAF-", name_ref, "-nonambiguous.png")
          note_type <- "non-ambiguous SNPs"
          note_correlation1 <- "Non-ambiguous SNPs allele frequencies correlate poorly with those in"
          note_correlation2 <- " - - warning: allele frequency of non-ambiguous SNPs correlates poorly with reference (r = "
          note_failure <- "insufficient non-missing allele frequencies to correlate non-ambiguous SNPs allele frequencies with reference"
        }
      }
      
      if(sum( !is.na(dat_FRQ) & !is.na(ref_FRQ) ) > 10L ) {
        cor_FRQ <- cor(dat_FRQ, ref_FRQ, use = "na.or.complete")
        if(cor_FRQ < threshold_r) {
          if(use_log) save_log(phaseL = 3L, checkL = paste("allele match -", name_ref, sep=" "), typeL = note_type, SNPL = relN, allSNPs = allN, actionL = "-", noteL = paste(note_correlation1, name_ref,"( r = ", round(cor_FRQ, digits = 3), ")"), fileL = paste(name_dir, name_file, sep = "/"))
          print(paste0(note_correlation2, round(cor_FRQ, digits = 3), ")"), quote = FALSE)
        }
        if(plot_FRQ & (cor_FRQ < threshold_r | !plot_if_threshold)) {
          
          png(plot_file, width = 720, height = 720, res = 144)
          if(plot_intensity) {
            intensity_plot(x = dat_FRQ, y = ref_FRQ, strata = HQ_subset,
                           xmax = 1, xmin = 0, ymax = 1, ymin = 0, verbose = FALSE,
                           main= paste(name_ref, "allele-frequency correlation"), sub = name_data, cex.sub = 1.3,
                           xlab="Reported allele frequency", ylab="Minor allele frequency")
            if(plot_HQ) legend(0, 0.94, c(" Low quality", "High quality", "Equal intensity"), pch = 19, col = c("red", "black", "turquoise3"))
          } else {
            plot(0.5, 0.5, xlim = c(0,1), ylim = c(0,1), col = "white",
                 main=paste(name_ref, "allele-frequency correlation"), sub = name_data, cex.sub = 1.3,
                 xlab="Reported allele frequency", ylab="Minor allele frequency")
            if(plot_HQ) points(dat_FRQ[!HQ_subset], ref_FRQ[!HQ_subset], pch = 20, col = "grey" , cex = 0.8)
            points(dat_FRQ[ HQ_subset], ref_FRQ[ HQ_subset], pch = 20, col = "black", cex = 0.8)
            if(plot_HQ) legend(0, 0.94, c(" Low quality", "High quality"), pch = 20, col = c("grey", "black"))
          }
          text(0, 0.98, paste("r =", round(cor_FRQ, digits = 3)), pos = 4, cex=1.2, col = ifelse(cor_FRQ < threshold_r, "red", "black") )
          if(!is.na(ambiguous)) mtext(note_type, line = 0.5, cex = 1.2)
          dev.off()
        }
      } else {
        if(use_log) { save_log(phaseL = 3L, checkL = paste("allele match -", name_ref), typeL = note_type, SNPL = relN, allSNPs = allN, actionL = "-", noteL = note_failure, fileL = paste(name_dir, name_file, sep = "/"))
        } else { print(paste(" - -", note_failure) , quote = FALSE) }
        cor_FRQ <- NA
      }
      return(cor_FRQ)
    }
    
    
    FRQcor <- plot_allele_cor(
      dat_FRQ = FRQ_all, ref_FRQ = dataset$MAF,
      HQ_subset = HQ_subset,
      plot_FRQ = plot_FRQ, plot_if_threshold = plot_if_threshold,
      threshold_r = threshold_r,
      plot_intensity = plot_intensity, plot_HQ = plot_HQ,
      name_dir = save_dir, name_file = save_name,
      name_ref = ref_name, name_data = dataname,
      ambiguous = NA, use_log = use_log, relN = log_SNPall, allN = log_SNPall)
    
    if(!is.na(FRQcor) & check_ambiguous & SNPn_ambiguous > 0L) {
      
      FRQcor_ambi <- plot_allele_cor(
        dat_FRQ = ifelse(ambiguous_list, FRQ_all, NA), ref_FRQ = dataset$MAF,
        HQ_subset = HQ_subset,
        plot_FRQ = plot_FRQ, plot_if_threshold = plot_if_threshold,
        threshold_r = threshold_r,
        plot_intensity = plot_intensity, plot_HQ = plot_HQ,
        name_dir = save_dir, name_file = save_name,
        name_ref = ref_name, name_data = dataname,
        ambiguous = TRUE, use_log = use_log, relN = SNPn_ambiguous, allN = log_SNPall)
      
      FRQcor_others <- plot_allele_cor(
        dat_FRQ = ifelse(ambiguous_list, NA, FRQ_all), ref_FRQ = dataset$MAF,
        HQ_subset = HQ_subset,
        plot_FRQ = plot_FRQ, plot_if_threshold = plot_if_threshold,
        threshold_r = threshold_r,
        plot_intensity = plot_intensity, plot_HQ = plot_HQ,
        name_dir = save_dir, name_file = save_name,
        name_ref = ref_name, name_data = dataname,
        ambiguous = FALSE, use_log = use_log, relN = log_SNPall - SNPn_ambiguous, allN = log_SNPall)
      
    } else {
      
      FRQcor_ambi <- NA
      FRQcor_others <- NA
    }
    
  } else {
    
    FRQcor <- NA
    FRQcor_ambi <- NA
    FRQcor_others <- NA
    SNPn_diffEAF <- NA
  }
  
  return(list(FRQ_cor = FRQcor, FRQ_cor_ambiguous = FRQcor_ambi, FRQ_cor_nonambi = FRQcor_others,
              n_SNPs = nrow(dataset), n_missing = SNPn_missing, n_missing_data = SNPn_missing_data, n_missing_ref = SNPn_missing_ref,
              n_negative_strand = SNPn_min, n_negative_switch = SNPn_min_SS, n_negative_mismatch = SNPn_min_MM,
              n_strandswitch = SNPn_switch, n_mismatch = SNPn_mismatch,
              n_flipped = SNPn_flip, n_ambiguous = SNPn_ambiguous, n_suspect = SNPn_suspect, n_diffEAF = SNPn_diffEAF,
              
              MARKER = if(return_SNPs | return_ref_values) dataset$MARKER else NULL,
              EFFECT_ALL = if(return_SNPs) dataset$EFFECT_ALL else NULL,
              OTHER_ALL = if(return_SNPs) dataset$OTHER_ALL else NULL,
              STRAND = if(check_strand & return_SNPs) dataset$STRAND else NULL,
              EFFECT = if(col_effect & return_SNPs) dataset$EFFECT else NULL,
              EFF_ALL_FREQ = if(check_FRQ & return_SNPs) dataset$EFF_ALL_FREQ else NULL,
              
              ref_MINOR = if(return_ref_values) dataset$MINOR else NULL,
              ref_MAJOR = if(return_ref_values) dataset$MAJOR else NULL,
              ref_MAF = if(check_FRQ & return_ref_values) dataset$MAF else NULL))
}
