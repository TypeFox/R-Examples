QC_series <-
function(data_files, datafile_tag, output_filenames,
                      dir_data = getwd(), dir_output = paste(dir_data, "QCGWASed", sep="/"), dir_references = dir_data,
                      header_translations, out_header = "standard",
                      allele_ref_std, allele_name_std, allele_ref_alt, allele_name_alt,
                      update_alt = FALSE, update_savename, update_as_rdata = FALSE, backup_alt = FALSE, 
                      plot_effectsizes = TRUE, lim_effectsizes = NULL,
                      plot_SE = TRUE, label_SE = TRUE, plot_SK = TRUE, label_SK = "outliers",
                      save_filtersettings = FALSE, ...) {
  # Step 0: housekeeping & checking arguments
  
  stopifnot(is.character(dir_data), length(dir_data) == 1L, !is.na(dir_data), nchar(dir_data) > 0L,
            is.character(dir_output), length(dir_output) == 1L, !is.na(dir_output), nchar(dir_output) > 0L,
            is.character(dir_references), length(dir_references) == 1L, !is.na(dir_references), nchar(dir_references) > 0L,
            is.logical(update_alt), length(update_alt) == 1L, !is.na(update_alt),
            is.logical(plot_effectsizes), length(plot_effectsizes) == 1L, !is.na(plot_effectsizes),
            is.logical(plot_SE), length(plot_SE) == 1L, !is.na(plot_SE),
            is.logical(label_SE), length(label_SE) == 1L, !is.na(label_SE),
            is.logical(plot_SK), length(plot_SK) == 1L, !is.na(plot_SK),
            is.logical(label_SK) | is.character(label_SK), length(label_SK) == 1L, !is.na(label_SK),
            is.logical(save_filtersettings), length(save_filtersettings) == 1L, !is.na(save_filtersettings))
  if(is.logical(label_SK)) { label_SK <- if(label_SK) "outliers" else "none" }
  stopifnot(label_SK %in% c("all", "outliers", "none"))
  if(plot_effectsizes & !is.null(lim_effectsizes)) { stopifnot(is.numeric(lim_effectsizes), length(lim_effectsizes) == 2L, !any(is.na(lim_effectsizes)), lim_effectsizes[2] > lim_effectsizes[1]) }
  if(!file.exists(dir_data) | !file.exists(dir_references)) { stop("Cannot find specified directory(s)") }
  if(missing(data_files)) {
    if(missing(datafile_tag)) {
      print("", quote=FALSE)
      print("Please specify either file names or file tag! No QC performed.", quote=FALSE)
      print("", quote = FALSE)
      return(invisible(FALSE))
    }
    if(!is.character(datafile_tag) | length(datafile_tag) != 1L) stop("'datafile_tag' must be a single character string")
    if(is.na(datafile_tag)) stop("'datafile_tag' cannot be NA")
    if(any(c("/","\\",":","*","?","\"","<",">","|","^","{","(","+","[","$") %in% unlist(strsplit(datafile_tag,"")))) {
      print("", quote = FALSE)
      print("Invalid character in file tag! No QC performed.", quote = FALSE)
      print("", quote = FALSE)
      return(invisible(FALSE))
    }
    data_files <- list.files(path = dir_data, pattern = datafile_tag)
  } else {
    if(!is.character(data_files) & length(data_files) > 0L) stop("'data_files' must be a vector of character strings")
    if(any(is.na(data_files))) stop("'data_files' cannot be NA")
    if(any(duplicated(data_files))) stop("The same filename has been specified twice in 'data_files'")
  }
  if(!all(sapply(paste(dir_data, data_files, sep = "/"), file.exists))) { stop("cannot find specified data_file(s) in the input directory") }
  if(missing(output_filenames)) {
    output_filenames <- rep(NA, length(data_files))
  } else {
    if(length(data_files) != length(output_filenames)) { stop("The specified number of output filenames does not match the number of input filenames") }
    stopifnot(is.character(output_filenames), !any(duplicated(output_filenames)),
              !any(is.na(output_filenames)), !any(output_filenames == ""))		
  }
  
  fileN <- length(data_files)
  
  print("", quote = FALSE)
  print("Starting Quality Control for:", quote = FALSE)
  print(paste(" - File", 1:fileN, ":", data_files), quote = FALSE)
  print("", quote = FALSE)
  
  # checking/loading header_translations
  if(length(header_translations) == 1L) {
    stopifnot(is.character(header_translations), nchar(header_translations) > 2L, file.exists(paste(dir_references, header_translations, sep = "/")))
    print(paste(" - Header translation - loading table:", header_translations), quote = FALSE)
    flush.console()
    header_translations <- read.table(paste(dir_references, header_translations, sep = "/"), stringsAsFactors = FALSE)
  } else {	print(" - Header translation - table specified by user", quote = FALSE) }
  if(!is.data.frame(header_translations) & !is.matrix(header_translations)) { stop("'header_translations' is not a table, matrix or dataframe") }
  if(ncol(header_translations) != 2L) { stop("'header_translations' does not have two columns") }
  if(any(duplicated(header_translations[ ,2]))) { stop("'header_translations' contains duplicated elements in column 2") }
  
  # checking/loading out_header translations
  if(length(out_header) == 1L) {
    stopifnot(is.character(out_header), nchar(out_header) > 2L)
    if(out_header %in% c("original", "standard", "GWAMA", "PLINK", "META", "GenABEL", "old")) {
      print(paste(" - Output header - using option:", out_header), quote = FALSE)
    } else {
      if(!file.exists(paste(dir_references, out_header, sep = "/"))) { stop("'out_header' is not a standard name nor a filename in dir_references") }
      print(paste(" - Output header - loading table:", out_header), quote = FALSE)
      flush.console()
      out_header <- read.table(paste(dir_references, out_header, sep = "/"), stringsAsFactors = FALSE)
    }
  } else {
    if(is.matrix(out_header) | is.data.frame(out_header)) { print(" - Header translation - table specified by user", quote = FALSE)			
    } else { stop("'out_header' is not a table, filename or standard name") }
  }
  if(is.matrix(out_header) | is.data.frame(out_header)) {
    if(ncol(out_header) != 2L) { stop("'out_header' does not have two columns") }
    if(any(is.na(out_header))) { stop("'out_header' contains missing values") }
    if(any(duplicated(out_header[ ,1]),duplicated(out_header[ ,2]))) { stop("'out_header' contains duplicated names") }
  }
  
  remove_extension <- function(filename) {
    ci <- nchar(filename) - 1L
    ct <- "x"
    
    while(ct != "." & ci != 1L){
      ct <- substr(filename, ci, ci)
      ci <- ci - 1L
    }
    return( if(ci==1L) filename else substr(filename, 1L, ci) )
  }
  
  # checking/loading standard allele reference
  if(missing(allele_ref_std)) {
    allele_ref_std <- NULL
    allele_name_std <- "Standard"
  } else {
    if(!is.null(allele_ref_std)) {
      if(is.character(allele_ref_std)) {
        stopifnot(length(allele_ref_std) == 1L, nchar(allele_ref_std) > 2L, file.exists(paste(dir_references, allele_ref_std, sep = "/")))
        print(paste(" - Standard allele reference - loading table:", allele_ref_std), quote = FALSE)
        flush.console()
        if(toupper(substr(allele_ref_std, nchar(allele_ref_std) - 5L, nchar(allele_ref_std))) == ".RDATA") {
          if(missing(allele_name_std)) allele_name_std <- substr(allele_ref_std, 1L, nchar(allele_ref_std) - 6L)
          refs_name <- allele_ref_std
          refs_check <- load(paste(dir_references, allele_ref_std, sep = "/"))
          if(!"allele_ref_std" %in% refs_check) stop(paste0("no object with name 'allele_ref_std' inside ", refs_name,
                                                            ". The standard allele ref must have the object-name 'alele_ref_std', otherwise it cannot load."))
          if(length(refs_check) > 1L) print(paste(" - - WARNING:", refs_name, "contains more than one object. This may affect the running of the QC."), quote = FALSE)
        } else {
          if(missing(allele_name_std)) allele_name_std <- remove_extension(allele_ref_std)
          allele_ref_std <- read.table(paste(dir_references, allele_ref_std, sep = "/"), header = TRUE, stringsAsFactors = FALSE)
        }				
      } else {
        print(" - Standard allele reference - specified by user", quote = FALSE)
        if(missing(allele_name_std)) allele_name_std <- "Standard"
      }
      if(!is.data.frame(allele_ref_std) & !is.matrix(allele_ref_std)) stop("'allele_ref_std' isn't a table, matrix or dataframe")
      if(!all(c("SNP", "MAJOR", "MINOR", "MAF") %in% colnames(allele_ref_std))) {
        stop("Cannot find 'SNP', 'MAJOR', 'MINOR' or 'MAF' columns in 'allele_ref_std' table") }
      stopifnot(is.character(allele_name_std), length(allele_name_std) == 1L, !is.na(allele_name_std), nchar(allele_name_std) > 0L)
    } else { allele_name_std <- "Std. reference" }
  }
   
  # checking/loading alternative allele reference
  if(missing(allele_ref_alt)) {
    allele_ref_alt <- NULL
    if(missing(allele_name_alt)) { allele_name_alt <- "Alternative reference" }
  } else {	
    if(!is.null(allele_ref_alt)) {
      if(is.character(allele_ref_alt)) {
        stopifnot(length(allele_ref_alt) == 1L, nchar(allele_ref_alt) > 2L, file.exists(paste(dir_references, allele_ref_alt, sep = "/")))				
        print(paste(" - alternative allele reference - loading table:", allele_ref_alt), quote = FALSE)
        flush.console()
        if(toupper(substr(allele_ref_alt, nchar(allele_ref_alt) - 5L, nchar(allele_ref_alt))) == ".RDATA") {
          if(missing(update_savename)) update_savename <- substr(allele_ref_alt, 1L, nchar(allele_ref_alt) - 6L)
          if(missing(allele_name_alt)) allele_name_alt <- substr(allele_ref_alt, 1L, nchar(allele_ref_alt) - 6L)
          refa_name <- allele_ref_alt
          refa_check<- load(paste(dir_references, allele_ref_alt, sep = "/"))
          if(!"allele_ref_alt" %in% refa_check) stop(paste0("no object with name 'allele_ref_alt' inside ", refa_name,
                                                            ". The alternative allele reference must have the object-name 'alele_ref_alt', otherwise it cannot load."))
          if(length(refa_check) > 1L) print(paste(" - - WARNING:", refa_name, "contains more than one object. This may affect the running of the QC."), quote = FALSE)
        } else {
          if(missing(update_savename)) update_savename <- remove_extension(allele_ref_alt)
          if(missing(allele_name_alt)) allele_name_alt <- remove_extension(allele_ref_alt)
          allele_ref_alt <- read.table(paste(dir_references, allele_ref_alt, sep = "/"),
                                       sep = "\t", comment.char = "", header = TRUE,
                                       colClasses = c("character","integer","integer", "character","character","numeric","character","character"))
        }				
      } else {
        print(" - alternative allele reference - specified by user", quote = FALSE)
        if(missing(allele_name_alt)) allele_name_alt <- "Alternative reference"
        if(update_alt & missing(update_savename)) stop("No filename specified for saving updated allele-reference")
      }
      if(!is.data.frame(allele_ref_alt) & !is.matrix(allele_ref_alt)) stop("'allele_ref_alt' isn't a table, matrix or dataframe")
      if(!all(colnames(allele_ref_alt) == c("SNP", "CHR", "POS", "MINOR", "MAJOR", "MAF", "SOURCE", "DATE_ADDED"))) {
        stop("columns of 'allele_ref_alt' table do not match the standard") }
    } else {
      if(missing(allele_name_alt)) allele_name_alt <- "Alternative reference"
    }
  }
  
  if(update_alt) {
    stopifnot(is.character(update_savename), length(update_savename) == 1L, !is.na(update_savename), nchar(update_savename) > 0L,
              is.logical(update_as_rdata), length(update_as_rdata) == 1L, !is.na(update_as_rdata), is.logical(backup_alt), length(backup_alt) == 1L, !is.na(backup_alt))
  }
  
  if(!file.exists(dir_output)) {
    if(!dir.create(path = dir_output, showWarnings = TRUE, recursive = TRUE)) stop("Unable to create output directory")
    if(!file.exists(dir_output)) stop("Output directory was created but misnamed. Check the dir_output argument for unconventional characters or trailing spaces.")
  }
  print(paste0(" - Results saved in ", dir_output), quote=FALSE)	
  
  # Creating tables to store the results	
  QC_checklist <- data.frame(files =
                               c(	"Filename", "Samplesize N",
                                  "Lambda", "> Lambda - genotyped", "> Lambda - imputed",
                                  "SNPs in input file", "> monomorphic", ">> identical alleles *", "> excluded chromosomes",
                                  "SNPs pre-QC",	"> unusable",
                                  "SNPs mid-QC",	"> allele mismatch",
                                  "SNPs post-QC", "> % QC-removed **", "> % invalid", "> % genotyped", "> % imputed",
                                  "Fixed HWE p-value", "Fixed call rate", "Fixed sample size", "Fixed imputation quality",
                                  "Effect size 25%", "Effect size mean", "Effect size median", "Effect size 75%", "SE median", 
                                  "Allele frequency correlation", "Ambiguous allele frequency correlation", "p-value correlation", "Visscher's statistic ***",
                                  "Negative-strand SNPs", "Standard columns missing", "Standard columns empty", "Non-standard columns not converted",
                                  "", "Allele matching", "Reference strand switch", "> double strand switch", "Other strand switch", "> double strand switch",
                                  "Ambiguous SNPs", "> inverted allele frequency", "SNPs with deviating allele frequency", "", 
                                  " * Does not include SNPs with allele frequency = 0 or 1", " ** The % QC-removed is the % of pre-QC SNPs", " *** Calculated over high-quality SNPs only"
                               ), matrix(data = "", nrow = 48, ncol = fileN, dimnames = list(NULL, data_files)), stringsAsFactors = FALSE)
  table_plot	<- data.frame(no = c(1:fileN), dataset = data_files, QC_succes = FALSE, SE = 0, N = 0, skewness = 0, kurtosis = 0, HQ_SNPs = 0, stringsAsFactors = FALSE)
  if(plot_effectsizes) { table_ES <- data.frame(matrix(data = 0, nrow = 1000, ncol = fileN)) }
  if(save_filtersettings) {
    table_filter <- data.frame(file = data_files,
                               filter_FRQ = numeric(length = fileN), filter_HWE = numeric(length = fileN),
                               filter_cal = numeric(length = fileN), filter_imp = numeric(length = fileN),
                               FRQ_NA = logical(length = fileN), HWE_NA = logical(length = fileN),
                               cal_NA = logical(length = fileN), imp_NA = logical(length = fileN),
                               ignore_impstatus = logical(length = fileN), stringsAsFactors = FALSE) }
  all_ref_updated <- FALSE
  
  # Running the actual QC
  for (fileI in 1:fileN) {
    resultsI <- QC_GWAS(
      filename = data_files[fileI],
      filename_output = output_filenames[fileI],
      
      dir_data = dir_data, dir_output = dir_output, dir_references = dir_references,
      header_translations = header_translations, out_header = out_header,
      
      allele_ref_std = allele_ref_std, allele_name_std = allele_name_std,
      allele_ref_alt = allele_ref_alt, allele_name_alt = allele_name_alt,
      update_alt = update_alt, update_savename = update_savename,
      backup_alt = backup_alt, update_as_rdata = update_as_rdata,
      return_HQ_effectsizes = plot_effectsizes,
      logI = fileI, logN = fileN, ...)
    
    if(resultsI$QC_successful) {
      QC_checklist[ , fileI + 1] <- c(
        resultsI$filename_input, resultsI$sample_size,
        round(c(resultsI$lambda, resultsI$lambda_geno, resultsI$lambda_imp), digits = 3),
        resultsI$SNP_N_input, resultsI$SNP_N_input_monomorphic, resultsI$SNP_N_input_monomorphic_identic_alleles, resultsI$SNP_N_input_chr,
        resultsI$SNP_N_preQC, resultsI$SNP_N_preQC_unusable,
        resultsI$SNP_N_midQC, resultsI$SNP_N_midQC_mismatch,
        resultsI$SNP_N_postQC, round(100*(resultsI$SNP_N_preQC_unusable + resultsI$SNP_N_midQC_mismatch) / resultsI$SNP_N_preQC, digits = 2), round(100*c(resultsI$SNP_N_postQC_invalid, resultsI$SNP_N_postQC_geno, resultsI$SNP_N_postQC_imp) / resultsI$SNP_N_postQC, digits = 2),
        resultsI$fixed_HWE, resultsI$fixed_callrate, resultsI$fixed_sampleN, resultsI$fixed_impQ,
        resultsI$effect_25, resultsI$effect_mean, resultsI$effect_median, resultsI$effect_75,	resultsI$SE_median,
        round(c(resultsI$all_MAF_std_r, resultsI$all_ambiguous_MAF_std_r), digits = 3), round(resultsI$pvalue_r, digits = 3), round(resultsI$visschers_stat_HQ, digits = 3),
        if(is.na(resultsI$SNP_N_preQC_min)) { "no data" } else { resultsI$SNP_N_preQC_min > 0L },
        paste(resultsI$columns_std_missing, collapse = " "), paste(resultsI$columns_std_empty, collapse = " "), paste(resultsI$columns_unidentified, collapse = " "),
        "", "",
        resultsI$SNP_N_midQC_strandswitch_std, resultsI$SNP_N_midQC_strandswitch_std_min, resultsI$SNP_N_midQC_strandswitch_alt, resultsI$SNP_N_midQC_strandswitch_alt_min,
        resultsI$SNP_N_midQC_ambiguous, resultsI$SNP_N_midQC_suspect,resultsI$SNP_N_midQC_diffEAF,
        "", "", "", "" )
      if(resultsI$all_ref_changed) {
        if(!all_ref_updated) {
          update_savename_full <- paste0(dir_references, "/", update_savename, if(update_as_rdata) ".RData" else ".txt")
          backup_alt <- FALSE # so that only a single back-up is made
          all_ref_updated <- TRUE
        }
        if(fileI != fileN) { # so that ref_alt is not reloaded if there is no further QC
          print("", quote = FALSE)
          print("", quote = FALSE)
          print(" - - reloading the (updated) alternative allele reference", quote = FALSE)
          flush.console()
          if(update_as_rdata){
            load(update_savename_full)
          } else {
            allele_ref_alt <- read.table(update_savename_full, header = TRUE, sep = "\t", comment.char = "",
                                         colClasses = c("character","integer","integer", "factor","factor","numeric","factor","factor"))
          }
        }
      }
      
      if(plot_effectsizes) { table_ES[ , fileI] <- if(resultsI$effectsize_return) resultsI$effectsizes_HQ else NA }
      table_plot[fileI, "QC_succes"] <- TRUE
      table_plot[fileI, "SE"] <- resultsI$SE_median_HQ
      table_plot[fileI, "N" ] <- resultsI$sample_size_HQ
      table_plot[fileI, "skewness"] <- resultsI$skewness_HQ
      table_plot[fileI, "kurtosis"] <- resultsI$kurtosis_HQ
      table_plot[fileI, "HQ_SNPs"]	<- resultsI$SNP_N_postQC_HQ
      
      if(save_filtersettings) {
        table_filter$file[fileI] <- resultsI$filename_output    
        if(is.null(resultsI$settings_filter_HQ_FRQ)) {
          table_filter$filter_FRQ[fileI] <- NA
          table_filter$FRQ_NA[fileI] <- FALSE
        } else {
          table_filter$filter_FRQ[fileI] <- resultsI$settings_filter_HQ_FRQ
          table_filter$FRQ_NA[fileI] <- resultsI$settings_filter_NA_FRQ      
        }
        if(is.null(resultsI$settings_filter_HQ_HWE)) {
          table_filter$filter_HWE[fileI] <- NA
          table_filter$HWE_NA[fileI] <- FALSE
        } else {
          table_filter$filter_HWE[fileI] <- resultsI$settings_filter_HQ_HWE
          table_filter$HWE_NA[fileI] <- resultsI$settings_filter_NA_HWE      
        }
        if(is.null(resultsI$settings_filter_HQ_cal)) {
          table_filter$filter_cal[fileI] <- NA
          table_filter$cal_NA[fileI] <- FALSE
        } else {
          table_filter$filter_cal[fileI] <- resultsI$settings_filter_HQ_cal
          table_filter$cal_NA[fileI] <- resultsI$settings_filter_NA_cal      
        }
        if(is.null(resultsI$settings_filter_HQ_imp)) {
          table_filter$filter_imp[fileI] <- NA
          table_filter$imp_NA[fileI] <- FALSE
        } else {
          table_filter$filter_imp[fileI] <- resultsI$settings_filter_HQ_imp
          table_filter$imp_NA[fileI] <- resultsI$settings_filter_NA_imp      
        }
        table_filter$ignore_impstatus[fileI] <- resultsI$settings_ignore_impstatus
      }
    } else {
      QC_checklist[1, fileI + 1L] <- resultsI$filename_input
      QC_checklist[2, fileI + 1L] <- "FAILED QC"
      if(plot_effectsizes) table_ES[ , fileI] <- NA
      table_plot[fileI, 4:7] <- NA
    }
    if(fileI != fileN) gc(verbose = FALSE) # frees the memory used inside QC_GWAS & read.table
  }
  
  if(any(table_plot$QC_succes)) {
    print("", quote = FALSE)
    print("Step 6: Between study checks", quote = FALSE)
    print(" - Saving checklist", quote = FALSE)    
    write.table(QC_checklist, paste(dir_output, "Checklist.txt", sep = "/"), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
    write.table(table_plot, paste(dir_output, "Checkgraph_legenda.txt", sep = "/"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
    if(save_filtersettings) { write.table(table_filter[table_plot$QC_succes, ], paste(dir_output, "Check_filtersettings.txt", sep = "/"), row.names=FALSE, quote=FALSE, sep="\t",
                                          append = file.exists(paste(dir_output, "Check_filtersettings.txt", sep = "/")), col.names=!file.exists(paste(dir_output, "Check_filtersettings.txt", sep = "/")), ) }
    
    if (nrow(table_plot) > 1L ) {
      temp_names <- gsub("_", ".", table_plot$dataset)
      temp_names <- gsub("-", ".", temp_names)
      posdot <- regexpr("\\.", temp_names)
      label_names <- substr(temp_names, 1, posdot - 1)
      while (any(duplicated(label_names))) {
        label_dupli <- label_names %in% label_names[duplicated(label_names)]
        posnew <- posdot
        posnew[label_dupli] <- posdot[label_dupli] + regexpr("\\.", substr(temp_names[label_dupli], posdot[label_dupli]+1, nchar(temp_names[label_dupli])) ) 
        if(all(posnew <= posdot)) break
        posdot[label_dupli] <- ifelse(posnew[label_dupli] < posdot[label_dupli], nchar(temp_names)[label_dupli], posnew[label_dupli])
        label_names <- substr(temp_names, 1, ifelse(posdot == nchar(temp_names), posdot, posdot - 1L))
      }
      
      if(any(nchar(label_names) > 15)) label_names <- table_plot$no      
      
      if(plot_SE) {
        print(" - Creating scatterplot Standard Error vs Sample Size", quote = FALSE)
        flush.console()
        plot_precision(SE = table_plot$SE, N = table_plot$N, labels = if(label_SE) label_names else NULL,
                       save_name = "Checkgraph_precision", save_dir = dir_output,
                       sub = "High-quality SNPs only!", col.sub = "red", cex.sub = 1) }
      if(plot_SK) {
        print(" - Creating Skewness vs Kurtosis plot", quote = FALSE)
        flush.console()
        plot_skewness(skewness = table_plot$skewness, kurtosis = table_plot$kurtosis, labels = label_names, plot_labels = label_SK,
                      save_name = "Checkgraph_skew&kurt", save_dir = dir_output,
                      sub = "High-quality SNPs only!", col.sub = "red", cex.sub = 1,
                      xlim = c(min(table_plot$skewness, na.rm = TRUE), max(table_plot$skewness, na.rm = TRUE) + 0.3 * (max(table_plot$skewness, na.rm = TRUE) - min(table_plot$skewness, na.rm = TRUE)))) }
      if(plot_effectsizes) {
        print(" - Creating combined effect-size boxplot", quote = FALSE)
        flush.console()
        plot_distribution(data_table = table_ES, names = paste0(label_names, "; N =", table_plot$N),
                          include = table_plot$QC_succes, plot_order = table_plot$N,
                          quantile_lines = TRUE, save_name = "Checkgraph_effect-size", save_dir = dir_output,
                          ylim = lim_effectsizes, las = 0,
                          main = "Effect-size distribution", sub = "High-quality SNPs only!", col.sub = "red", cex.sub = 1)
  } } }
  
  if(plot_effectsizes) rm(table_ES)
  if(save_filtersettings) rm(table_filter)
  rm(header_translations, out_header, allele_ref_std, allele_ref_alt, QC_checklist, table_plot, resultsI)
  gc(verbose = FALSE)
  print("", quote = FALSE)
  print("", quote = FALSE)
  return(invisible(all_ref_updated))	
}
