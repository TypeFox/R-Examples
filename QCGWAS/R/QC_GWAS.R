QC_GWAS <-
function(filename, filename_output = paste0("QC_", filename),
                    dir_data = getwd(), dir_output = paste(dir_data, "QCGWASed", sep="/"), dir_references = dir_data,
                    header_translations,
                    column_separators = c("\t", " ", "", ",", ";"),
                    nrows = -1, nrows_test = 1000, header = TRUE, comment.char = "",
                    na.strings = c("NA", "nan", "NaN", "."), 
                    imputed_T = c("1", "TRUE", "T"), imputed_F = c("0", "FALSE", "F"), imputed_NA = c(NA, "-"),

                    save_final_dataset = TRUE, gzip_final_dataset = TRUE, order_columns = FALSE,
                    spreadsheet_friendly_log = FALSE,
                    out_header = "standard", out_quote = FALSE, out_sep = "\t",
                    out_eol = "\n", out_na = "NA", out_dec = ".", out_qmethod = "escape",
                    out_rownames = FALSE, out_colnames = TRUE,
                    return_HQ_effectsizes = FALSE,
                    
                    remove_X = FALSE, remove_Y = FALSE, remove_XY = remove_Y, remove_M = FALSE,
                    calculate_missing_p = FALSE,
                    
                    make_plots = TRUE, only_plot_if_threshold = TRUE,
                    threshold_allele_freq_correlation = 0.95, threshold_p_correlation = 0.99,
                    plot_intensity = FALSE,
                    plot_histograms = make_plots,  plot_QQ = make_plots, plot_QQ_bands = TRUE,
                    plot_Manhattan = make_plots, plot_cutoff_p = 0.05,
                    
                    allele_ref_std, allele_name_std,
                    allele_ref_alt, allele_name_alt,
                    update_alt = FALSE, update_savename, update_as_rdata = FALSE, backup_alt = FALSE,
                    remove_mismatches = TRUE,
                    remove_mismatches_std = remove_mismatches, remove_mismatches_alt = remove_mismatches,
                    threshold_diffEAF = 0.15, remove_diffEAF = FALSE,
                    remove_diffEAF_std = remove_diffEAF, remove_diffEAF_alt = remove_diffEAF,
                    check_ambiguous_alleles = FALSE, 
                    
                    use_threshold = 0.1,
                    useFRQ_threshold = use_threshold, useHWE_threshold = use_threshold, useCal_threshold = use_threshold, useImp_threshold = use_threshold, useMan_threshold = use_threshold,
                    HQfilter_FRQ = 0.01, HQfilter_HWE = 10^-6, HQfilter_cal = 0.95, HQfilter_imp = 0.3,
                    QQfilter_FRQ = c(NA, 0.01, 0.05), QQfilter_HWE = c(NA, 10^-6, 10^-4),
                    QQfilter_cal = c(NA, 0.95, 0.99), QQfilter_imp = c(NA, 0.3, 0.5, 0.8),
                    NAfilter = TRUE, NAfilter_FRQ = NAfilter,
                    NAfilter_HWE = NAfilter, NAfilter_cal = NAfilter, NAfilter_imp = NAfilter,
                    ignore_impstatus = FALSE,
                    minimal_impQ_value = -0.5, maximal_impQ_value = 1.5,
                    logI = 1L, logN = 1L, ...) {
  
  ### PHASE 0: checking the input parameters for errors
  # Filename_output is checked later, to avoid calling it before the extensions
  #	are removed from "filename". If filename_output is invalid, "QC_[filename]"
  #	will be used instead.
  
  start_time <- date()
  
  stopifnot(is.character(filename), length(filename) == 1L, !is.na(filename),
            is.character(dir_data), length(dir_data) == 1L,
            is.character(dir_output), length(dir_output) == 1L,
            is.character(dir_references), length(dir_references) == 1L)
  if(nchar(filename) < 3L | nchar(dir_data) == 0L | nchar(dir_output) == 0L | nchar(dir_references) == 0L) {
    stop("File / directory names cannot be empty character-strings") } 
  if(!file.exists(dir_data) | !file.exists(dir_references)) stop("Cannot find specified directory(s)")
  if(!file.exists(paste(dir_data, filename, sep = "/"))) stop("Cannot find the data file in input directory")
  
  if(length(header_translations) == 1L) {
    stopifnot(is.character(header_translations), nchar(header_translations) > 2L, file.exists(paste(dir_references, header_translations, sep = "/")))
    settings_header_input <- header_translations
    header_translations <- read.table(paste(dir_references, header_translations, sep = "/"), stringsAsFactors = FALSE)
  } else { settings_header_input <- "table" }
  if(!is.data.frame(header_translations) & !is.matrix(header_translations)) stop("'header_translations' is not a table, matrix or dataframe")
  if(ncol(header_translations) != 2L) stop("'header_translations' does not have two columns")
  if(any(duplicated(header_translations[ ,2]))) stop("'header_translations' contains duplicated elements in column 2")
  
  stopifnot(is.logical(header), length(header) == 1L)
  if(is.na(header))	stop("'header' cannot be NA")
  stopifnot(is.numeric(nrows), length(nrows) == 1L)
  if(is.na(nrows) | nrows == 0)	stop("'nrows' cannot be missing or zero")
  stopifnot(is.numeric(nrows_test), length(nrows_test) == 1L)
  if(is.na(nrows_test) | nrows_test == 0) stop("'nrows_test' cannot be missing or zero") 
  stopifnot(is.character(comment.char), length(comment.char) == 1L)
  if(nchar(comment.char) != 1L & nchar(comment.char) != 0L) stop("'comment.char' must be a single character or empty string")
  
  stopifnot(is.vector(na.strings), is.character(na.strings), is.vector(imputed_T), is.vector(imputed_F), is.vector(imputed_NA),
            is.vector(column_separators), is.character(column_separators), is.logical(ignore_impstatus), length(ignore_impstatus) == 1L, !is.na(ignore_impstatus))
  if(any(duplicated(c(imputed_NA, imputed_T, imputed_F)))) stop("duplicate strings in the 'imputed' arguments")
  
  stopifnot(is.logical(save_final_dataset), length(save_final_dataset) == 1L, is.logical(order_columns), length(order_columns) == 1L,
            is.logical(spreadsheet_friendly_log), length(spreadsheet_friendly_log) == 1L)
  if(is.na(save_final_dataset)) stop("'save_final_dataset' cannot be NA")
  if(is.na(order_columns)) stop("'order_columns' cannot be NA")
  if(is.na(spreadsheet_friendly_log)) stop("'spreadsheet_friendly_log' cannot be NA")
  
  if(save_final_dataset) {
    stopifnot(is.logical(out_quote), length(out_quote) == 1L,
              is.character(out_sep), length(out_sep) == 1L,
              is.character(out_eol), length(out_eol) == 1L,
              is.character(out_na),	length(out_na) == 1L,
              is.character(out_dec), length(out_dec) == 1L,
              is.character(out_qmethod), length(out_qmethod) == 1L,
              is.logical(gzip_final_dataset), length(gzip_final_dataset) == 1L)
    if(is.na(out_quote)) stop("'out_quote' cannot be NA")
    if(is.na(gzip_final_dataset)) stop("'gzip_final_dataset' cannot be NA")
    if(nchar(out_sep) == 0L) stop("'out_sep' cannot be an empty character-string")
    if(nchar(out_eol) == 0L) stop("'out_eol' cannot be an empty character-string")
    if(nchar(out_dec) != 1L) stop("'out_dec' must be of length 1")
    if(!out_qmethod %in% c("escape", "double", "e", "d")) stop("'out_qmethod' must be either 'escape' or 'double'")
    if(length(out_header) == 1L) {
      stopifnot(is.character(out_header), nchar(out_header) > 2L)
      settings_header_output <- out_header
      if(out_header %in% c("original", "standard", "GWAMA", "PLINK", "META", "GenABEL", "old")) {
        if(settings_header_output != "standard" & settings_header_output != "original") {
          if(settings_header_output == "GenABEL"){
            out_header <- data.frame(
              GenABEL = c("name", "chromosome", "position", "strand", "allele1", "allele2", "effallelefreq", "n", "beta", "sebeta", "p", "pexhwe", "call"),
              QC = c("MARKER", "CHR", "POSITION", "STRAND", "EFFECT_ALL", "OTHER_ALL", "EFF_ALL_FREQ", "N_TOTAL", "EFFECT", "STDERR", "PVALUE", "HWE_PVAL", "CALLRATE"),
              stringsAsFactors = FALSE)
          } else {
            out_header <- data.frame(
              GWAMA = c("MARKER", "CHR", "POSITION", "EA", "NEA", "STRAND", "BETA", "SE", "P", "EAF", "N", "IMPUTED", "IMP_QUALITY"),
              PLINK = c("SNP",		"CHR", "BP",			 "A1", "A2",	"STRAND", "BETA", "SE", "P", "EFF_ALL_FREQ", "N", "IMPUTED", "IMP_QUALITY"),
              META = c( "rsid",	 "chr", "pos",			"allele_B", "allele_A", "strand", "beta", "se", "P_value", "EFF_ALL_FREQ", "N", "imputed", "info"),
              old =c("MARKER", "CHR", "POSITION", "ALLELE1",		"ALLELE2",	 "STRAND", "EFFECT", "STDERR", "PVALUE", "FREQLABEL",		"N_TOTAL", "IMPUTED", "IMP_QUALITY"),
              QC = c("MARKER", "CHR", "POSITION", "EFFECT_ALL", "OTHER_ALL", "STRAND", "EFFECT", "STDERR", "PVALUE", "EFF_ALL_FREQ", "N_TOTAL", "IMPUTED", "IMP_QUALITY"),
              stringsAsFactors = FALSE)
            out_header <- out_header[ ,c(settings_header_output, "QC")]
          }
        }
      } else {
        if(!file.exists(paste(dir_references, out_header, sep = "/"))) stop("'out_header' is not a standard name nor a filename in dir_references")
        out_header <- read.table(paste(dir_references, out_header, sep = "/"), stringsAsFactors = FALSE)
      }
    } else {
      if(!is.matrix(out_header) & !is.data.frame(out_header)) stop("'out_header' is not a table, filename or standard name")
      settings_header_output <- "table"
    }
    if(is.matrix(out_header) | is.data.frame(out_header)) {
      if(ncol(out_header) != 2L) { stop("'out_header' does not have 2 columns") }
      if(any(is.na(out_header))) { stop("'out_header' contains missing values") }
      if(any(duplicated(out_header[ ,1]),duplicated(out_header[ ,2]))) { stop("'out_header' contains duplicated names") }
    }
  } else { settings_header_output <- NA }
  
  stopifnot(is.logical(remove_X), length(remove_X) == 1L,
            is.logical(remove_Y), length(remove_Y) == 1L,
            is.logical(remove_XY), length(remove_XY) == 1L,
            is.logical(remove_M), length(remove_M) == 1L)
  if(is.na(remove_X)) { stop( "'remove_X' cannot be NA") }
  if(is.na(remove_Y)) { stop( "'remove_Y' cannot be NA") }
  if(is.na(remove_XY)){ stop("'remove_XY' cannot be NA") }
  if(is.na(remove_M)) { stop( "'remove_M' cannot be NA") }
  
  stopifnot(is.numeric(useFRQ_threshold), length(useFRQ_threshold) == 1L,
            is.numeric(useHWE_threshold), length(useHWE_threshold) == 1L,
            is.numeric(useCal_threshold), length(useCal_threshold) == 1L,
            is.numeric(useImp_threshold), length(useImp_threshold) == 1L,
            is.numeric(useMan_threshold), length(useMan_threshold) == 1L)
  if(any(is.na(c(useFRQ_threshold, useHWE_threshold, useCal_threshold, useImp_threshold, useMan_threshold))) | useFRQ_threshold <= 0 | useHWE_threshold <= 0 | useCal_threshold <= 0 | useImp_threshold <= 0 | useMan_threshold <= 0) {
    stop("'use_threshold' values cannot be negative, missing or zero") } 
  
  stopifnot(is.logical(plot_histograms), length(plot_histograms) == 1L,
            is.numeric(threshold_p_correlation), length(threshold_p_correlation) == 1L,
            is.logical(plot_QQ), length(plot_QQ) == 1L,
            is.logical(plot_QQ_bands), length(plot_QQ_bands) == 1L,
            is.logical(plot_Manhattan), length(plot_Manhattan) == 1L,
            is.numeric(plot_cutoff_p), length(plot_cutoff_p) == 1L)
  if(is.na(plot_histograms)) { stop("'plot_histograms' cannot be NA") }
  if(is.na(plot_QQ)) { stop("'plot_QQ' cannot be NA") }
  if(is.na(plot_QQ_bands)) { stop("'plot_QQ_bands' cannot be NA") }
  if(is.na(plot_Manhattan)) { stop("'plot_Manhattan' cannot be NA") }
  
  stopifnot(is.logical(NAfilter_FRQ), length(NAfilter_FRQ) == 1L,
            is.logical(NAfilter_HWE), length(NAfilter_HWE) == 1L,
            is.logical(NAfilter_cal), length(NAfilter_cal) == 1L,
            is.logical(NAfilter_imp), length(NAfilter_imp) == 1L)	
  if(is.na(NAfilter_FRQ) | is.na(NAfilter_HWE) | is.na(NAfilter_cal) | is.na(NAfilter_imp)) {
    stop("'NAfilter' values cannot be NA") }
  
  if(!is.null(HQfilter_FRQ)) {
    if(length(HQfilter_FRQ) > 1L) stop("'HQfilter_FRQ' isn't a single value")
    if(!is.numeric(HQfilter_FRQ) & !is.na(HQfilter_FRQ)) stop("'HQfilter_FRQ' isn't a numerical or NA value") }
  if(!is.null(HQfilter_HWE)) {
    if(length(HQfilter_HWE) > 1L) stop("'HQfilter_HWE' isn't a single value")
    if(!is.numeric(HQfilter_HWE) & !is.na(HQfilter_HWE)) stop("'HQfilter_HWE' isn't a numerical or NA value") }
  if(!is.null(HQfilter_cal)) {
    if(length(HQfilter_cal) > 1L) stop("'HQfilter_cal' isn't a single value")
    if(!is.numeric(HQfilter_cal) & !is.na(HQfilter_cal)) stop("'HQfilter_cal' isn't a numerical or NA value") }
  if(!is.null(HQfilter_imp)) {
    if(length(HQfilter_imp) > 1L) stop("'HQfilter_imp' isn't a single value")
    if(!is.numeric(HQfilter_imp) & !is.na(HQfilter_imp)) stop("'HQfilter_imp' isn't a numerical or NA value") }
  
  if(!is.numeric(minimal_impQ_value) | length(minimal_impQ_value) != 1L) stop("'minimal_impQ_value' values isn't a single numerical value")
  if(!is.numeric(maximal_impQ_value) | length(maximal_impQ_value) != 1L) stop("'maximal_impQ_value' values isn't a single numerical value")
  
  if(!is.null(QQfilter_FRQ)) { 
    if(!is.vector( QQfilter_FRQ)) stop("'QQfilter_FRQ' is not a numeric vector")
    if(!is.numeric(QQfilter_FRQ) & !all(is.na(QQfilter_FRQ))) stop("'QQfilter_FRQ' is not a numeric vector") }
  if(!is.null(QQfilter_HWE)) {
    if(!is.vector( QQfilter_HWE)) stop("'QQfilter_HWE' is not a numeric vector")
    if(!is.numeric(QQfilter_HWE) & !all(is.na(QQfilter_HWE))) stop("'QQfilter_HWE' is not a numeric vector") }
  if(!is.null(QQfilter_cal)) {
    if(!is.vector( QQfilter_cal)) stop("'QQfilter_cal' is not a numeric vector")
    if(!is.numeric(QQfilter_cal) & !all(is.na(QQfilter_cal))) stop("'QQfilter_cal' is not a numeric vector") }
  if(!is.null(QQfilter_imp)) {
    if(!is.vector( QQfilter_imp)) stop("'QQfilter_imp' is not a numeric vector")
    if(!is.numeric(QQfilter_imp) & !all(is.na(QQfilter_imp))) stop("'QQfilter_imp' is not a numeric vector") }
  
  stopifnot(is.logical(calculate_missing_p), length(calculate_missing_p) == 1L,
            is.numeric(logI), length(logI) == 1L, is.numeric(logN), length(logN) == 1L)
  if(is.na(calculate_missing_p)) { stop("'calculate_missing_p' cannot be NA") }
  
  use_allele_std <- if(missing(allele_ref_std)) FALSE else !is.null(allele_ref_std)
  use_allele_alt <- if(missing(allele_ref_alt)) FALSE else !is.null(allele_ref_alt)
  
  remove_extension <- function(filename) {
    ci <- nchar(filename) - 1L
    ct <- "x"
    while(ct != "." & ci != 1L){
      ct <- substr(filename, ci, ci)
      ci <- ci - 1L
    }
    return( if(ci==1L) filename else substr(filename, 1L, ci) )
  }
  
  if(!is.logical(update_alt) | length(update_alt) != 1L) { stop("'update_alt' is not a single logical value") }
  if(is.na(update_alt)) { stop("'update_alt' cannot be NA") }
  
  # arguments passed on to match_alleles are now tested regardless whether
  # the function is called, since they are used in the settings-tables of
  # the log file.
  stopifnot(is.logical(check_ambiguous_alleles), length(check_ambiguous_alleles) == 1L,
            is.numeric(threshold_allele_freq_correlation), length(threshold_allele_freq_correlation) == 1L,
            is.logical(remove_mismatches_std), length(remove_mismatches_std) == 1L,
            is.logical(remove_mismatches_alt), length(remove_mismatches_alt) == 1L,
            is.logical(remove_diffEAF_std), length(remove_diffEAF_std) == 1L,
            is.logical(remove_diffEAF_alt), length(remove_diffEAF_alt) == 1L,
            is.numeric(threshold_diffEAF), length(threshold_diffEAF) == 1L,
            is.logical(plot_intensity), length(plot_intensity) == 1L)
  if(is.na(check_ambiguous_alleles)) stop("'check_ambiguous_alleles' cannot be NA")
  if(is.na(remove_mismatches_std)) stop("'remove_mismatches_std' cannot be NA")
  if(is.na(remove_mismatches_alt)) stop("'remove_mismatches_alt' cannot be NA")
  if(is.na(remove_diffEAF_std)) stop("'remove_diffEAF_std' cannot be NA")
  if(is.na(remove_diffEAF_alt)) stop("'remove_diffEAF_alt' cannot be NA")
  if(is.na(plot_intensity)) stop("'plot_intensity' cannot be NA")
  
  if(use_allele_std) {
    if(is.character(allele_ref_std)) {
      stopifnot(length(allele_ref_std) == 1L, nchar(allele_ref_std) > 2L, file.exists(paste(dir_references, allele_ref_std, sep = "/")))
      print(paste("Loading standard allele reference from file:", allele_ref_std), quote = FALSE)
      flush.console()
      settings_allele_std <- allele_ref_std
      if(toupper(substr(allele_ref_std, nchar(allele_ref_std) - 5L, nchar(allele_ref_std))) == ".RDATA") {
        if(missing(allele_name_std)) allele_name_std <- substr(allele_ref_std, 1L, nchar(allele_ref_std) - 6L)
        refs_name <- allele_ref_std
        refs_check<- load(paste(dir_references, allele_ref_std, sep = "/"))
        if(!"allele_ref_std" %in% refs_check) stop(paste0("no object with name 'allele_ref_std' inside ", refs_name,
                                                          ". The standard allele reference must have the object-name 'alele_ref_std', otherwise it cannot load."))
        if(length(refs_check) > 1L) print(paste(" - - WARNING:", refs_name, "contains more than one object. This may affect the running of the QC."), quote = FALSE)
      } else {
        if(missing(allele_name_std)) allele_name_std <- remove_extension(allele_ref_std)
        allele_ref_std <- read.table(paste(dir_references, allele_ref_std, sep = "/"), header = TRUE, stringsAsFactors = FALSE)
      }
    } else {
      settings_allele_std <- "table"
      if(missing(allele_name_std)) allele_name_std <- "Std. reference"
    }
    if(!is.data.frame(allele_ref_std) & !is.matrix(allele_ref_std)) stop("'allele_ref_std' isn't a table, matrix or dataframe")
    if(!all(c("SNP", "MAJOR", "MINOR", "MAF") %in% colnames(allele_ref_std))) {
      stop("Cannot find 'SNP', 'MAJOR', 'MINOR' or 'MAF' columns in 'allele_ref_std' table") }
    stopifnot(is.character(allele_name_std), length(allele_name_std) == 1L, !is.na(allele_name_std), nchar(allele_name_std) > 0L)
  } else {
    settings_allele_std <- NA
    allele_name_std <- "Std. reference"
    remove_mismatches_std <- FALSE
    remove_diffEAF_std <- FALSE
  }
  
  if(use_allele_alt) {
    if(is.character(allele_ref_alt)) {
      stopifnot(length(allele_ref_alt) == 1L, nchar(allele_ref_alt) > 2L, file.exists(paste(dir_references, allele_ref_alt, sep = "/")))
      print(paste("Loading alternative allele reference from file:", allele_ref_alt), quote = FALSE)
      flush.console()
      settings_allele_alt <- allele_ref_alt
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
      settings_allele_alt <- "table"
      if(missing(allele_name_alt)) allele_name_alt <- "alternative reference"
      if(update_alt & missing(update_savename)) stop("No filename specified for saving updated allele-reference")
    }
    if(!is.data.frame(allele_ref_alt) & !is.matrix(allele_ref_alt)) stop("'allele_ref_alt' isn't a table, matrix or dataframe")
    if(!all(colnames(allele_ref_alt) == c("SNP", "CHR", "POS", "MINOR", "MAJOR", "MAF", "SOURCE", "DATE_ADDED"))) stop("columns of 'allele_ref_alt' table do not match the standard")
    if(!is.character(allele_name_alt) | length(allele_name_alt) != 1L) stop("'allele_name_alt' is not a character-string")
    if(nchar(allele_name_alt) == 0L) stop("'allele_name_alt' cannot be an empty character-string")
  } else {
    if(update_alt & missing(update_savename)) stop("Cannot update alternative reference: update_savename specified")
    settings_allele_alt <- NA
    allele_name_alt <- "alternative reference"
    remove_mismatches_alt <- FALSE
    remove_diffEAF_alt <- FALSE
  }
  
  if(update_alt) {
    stopifnot(is.character(update_savename), length(update_savename) == 1L, !is.na(update_savename), nchar(update_savename) > 0L,
              is.logical(update_as_rdata), length(update_as_rdata) == 1L, !is.na(update_as_rdata),
              is.logical(backup_alt), length(backup_alt) == 1L, !is.na(backup_alt))
  }
  
  if(!file.exists(dir_output)) {
    if(!dir.create(path = dir_output, showWarnings = TRUE, recursive = TRUE)) stop("Unable to create output directory")
    if(!file.exists(dir_output)) stop("Output directory was created but misnamed. Check the dir_output argument for unconventional characters or trailing spaces.")
  }
  
  ### PHASE 1: preparing the data
  # Step 1a: loading dataset. Two control variables are created here.
  #	EmergencyExit is a failsafe: when the script encounters a problem
  #	that cannot be fixed or ignored (corrupt data or when a QC step
  #	removes all SNPs in the dataset), EmergencyExit is set to TRUE. The
  #	script tests EmergencyExit in several places. If TRUE, the rest of
  #	the QC is skipped. allele_ref_changed keeps track of whether the
  #	alternative allele reference is changed.
  
  emergency_return <- function(name_in, name_out, LI, LN) {
    print("", quote = FALSE)
    print(paste("QC check aborted for", name_in,"( file", LI, "out of", LN, ")"), quote = FALSE)
    return(list(QC_successful = FALSE, filename_input = name_in, filename_output = name_out, all_ref_changed = FALSE, effectsize_return = FALSE))
  }
  
  EmergencyExit <- FALSE
  allele_ref_changed <- FALSE
  filename_input <- filename
  
  print("", quote = FALSE)
  print("", quote = FALSE)
  print(paste("Starting analysis of", filename_input,"( file", logI, "out of", logN, ")"), quote = FALSE)
  print(	"Step 1: loading dataset", quote = FALSE)
  flush.console()
  
  loadI <- load_test(filename_input, dir_data, column_separators = column_separators, test_nrows = nrows_test,
                     header = header, na.strings = na.strings, comment.char = comment.char, stringsAsFactors = FALSE, ...)
  if(loadI$success) {
    if(loadI$type == "zip" | loadI$type == ".gz") {
      if(loadI$type == "zip") {
        filename <- substr(filename, 1L, nchar(filename) - 4L)
        dataI <- read.table(unz(paste(dir_data, filename_input, sep = "/"), filename), sep = loadI$sep,
                            header = header, nrows = nrows, na.strings = na.strings, comment.char = comment.char, stringsAsFactors = FALSE, ...)
        close(unz(paste(dir_data, filename_input, sep = "/"), filename))
      } else {
        filename <- substr(filename, 1L, nchar(filename) - 3L)
        dataI <- read.table(gzfile(paste(dir_data, filename_input, sep = "/")), sep = loadI$sep,
                            header = header, nrows = nrows, na.strings = na.strings, comment.char = comment.char, stringsAsFactors = FALSE, ...)
        close(gzfile(paste(dir_data, filename_input, sep = "/")))
      }
    } else {
      dataI <- read.table(paste(dir_data, filename_input, sep = "/"), sep = loadI$sep,
                          header = header, nrows = nrows, na.strings = na.strings, comment.char = comment.char, stringsAsFactors = FALSE, ...)
    }
  } else { EmergencyExit <- TRUE }
  
  # Filename_output is checked here to avoid calling it before the extension has been removed from filename_input
  filename <- remove_extension(filename)
  filename_error <- TRUE
  
  if(missing(filename_output)) {		
    filename <- paste("QC", filename, sep = "_")
    filename_error <- FALSE
  } else {
    if(length(filename_output) == 1L) {
      if(is.na(filename_output)) {
        filename <- paste("QC", filename, sep = "_")
        filename_error <- FALSE
      } else {
        if(is.character(filename_output)) { if(nchar(filename_output) > 0L) {
          filename <- filename_output
          filename_error <- FALSE
  } } } } }
  
  if(filename_error) { filename <- paste("QC", filename, sep = "_") }
  filename_dir <- paste(dir_output, filename, sep = "/")
  
  write.table(t(c("Step", "Check", "Type", "Affected SNPs", "%", "Action", "Notes")), paste0(filename_dir, "_log.txt"), append=FALSE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
  
  if(filename_error) {
    save_log(1L, "loading dataset", "output filename", 0L, 1L, "-", "Cannot use specified output filename - reverted to default", filename_dir)
    print(paste0(" - - warning: unable to use specified output filename - reverted to default: ", filename, ".txt"), quote = FALSE)
  }
  if(EmergencyExit) { # 1st phase-1 EmergencyExit: failure to load file
    print("CRITICAL ERROR: unable to load dataset ", quote=FALSE)
    print(paste(" - error type:", loadI$error), quote=FALSE)
    save_log(1L, "loading data", loadI$error, 0L, 1L, "QC aborted", "Cannot read file format", filename_dir)
    return(emergency_return(name_in = filename_input, name_out = filename, LI = logI, LN = logN))
  }
  
  # Step 1b: checking and translating columnnames
  
  SNPn_input	<- nrow(dataI)
  header_orig <- colnames(dataI)
  header_info <- translate_header(header = header_orig, alternative = header_translations)
  
  # 2st phase-1 EmergencyExit: duplicate or missing crucial headers
  if(any(duplicated(header_info$header_h))) {
    print("CRITICAL ERROR: dataset contains duplicate data columns",quote=FALSE)
    print(paste(" - duplicate columns:", paste(header_info$header_h[duplicated(header_info$header_h)], collapse = ", ") ), quote = FALSE)
    save_log(1L, "loading file", "duplicate columns", SNPn_input, SNPn_input, "QC aborted", paste("Dataset contains duplicate columns for", paste(header_info$header_h[duplicated(header_info$header_h)], collapse = ", ") ), fileL = filename_dir)
    return(emergency_return(name_in = filename_input, name_out = filename, LI = logI, LN = logN))
  }
  
  if(header_info$missing_N > 0L) {
    if(any(c("MARKER", "EFFECT_ALL", "OTHER_ALL", "EFFECT", "STDERR") %in% header_info$missing_h)) {
      print(paste("CRITICAL ERROR: column(s)", paste(header_info$missing_h, collapse = ", "), "not found in dataset"), quote = FALSE)
      save_log(1L, "loading file", "missing column(s)", SNPn_input, SNPn_input, "QC aborted", paste("Column(s)", paste(header_info$missing_h, collapse = ", "), "not found in dataset."), filename_dir)
      if(header_info$unknown_N > 0L) { save_log(1L, "loading file", "unidentified column", SNPn_input, SNPn_input, "-", paste("Column(s)", paste(header_info$unknown_h, collapse = ", "), "were not found in the translation table."), filename_dir) }
      return(emergency_return(name_in = filename_input, name_out = filename, LI = logI, LN = logN))
    }
    print(paste(" - - warning: column(s)", paste(header_info$missing_h, collapse = ", "), "not found"), quote = FALSE)
    save_log(1L, "loading file", "missing column", SNPn_input, SNPn_input, "dummy column created", paste("Column(s)", paste(header_info$missing_h, collapse = ", "), "not found in dataset."), filename_dir)
    if(header_info$unknown_N > 0L) { save_log(1L, "loading file", "unidentified column", SNPn_input, SNPn_input, "-", paste("Column(s)", paste(header_info$unknown_h, collapse = ", "), "were not found in the translation table."), filename_dir) }
    dataI <- cbind(dataI, matrix(data = NA, nrow = SNPn_input, ncol = header_info$missing_N,
                                 dimnames = list(NULL, header_info$missing_h))) 
    colnames(dataI)[1:header_info$header_N] <- header_info$header_h
  } else { colnames(dataI) <- header_info$header_h }

  if(order_columns) { dataI <- dataI[ , c("MARKER","CHR","POSITION","EFFECT_ALL","OTHER_ALL","STRAND","EFFECT","STDERR","PVALUE","EFF_ALL_FREQ","HWE_PVAL","CALLRATE","N_TOTAL","IMPUTED", "USED_FOR_IMP", "IMP_QUALITY", header_info$unknown_h) ] }
  
  
  ### PHASE 2: checking data-integrity
  # Checking for invalid data entries and selecting SNPs for exclusion
  # Step 2a is removing the monomorphic SNPs (defined as SNPs with a
  #	missing or invalid non-effect allele, allele freq = 1 or 0, or identical
  #	alleles). This step is caried out first so that the monomorphic
  #	SNPs won't generate unnecesarry logs in phase 2c-d. Y & M-
  #	chromosome SNPs are removed here as well.
  # Step 2b is making a list of poorly-imputed SNPs (imputation quality > 0.3). Missing
  #	crucial variables for these SNPs won't be counted in the log
  #	file (but they will be shown in the stat file). This step also
  #	compiles lists of which SNPs are imputed and genotyped
  # Step 2c is checking the remaining columns for missing and invalid
  #	entries.
  # Step 2d is the removal of SNPs from the QC. SNPs are removed if they
  #	are unusuable (missing or invalid crucial variables). Crucial
  #	variables are: marker name, effect allele, effect size and standard
  #	error (bad non-effect alleles have already been removed in phase 2a).
  #	Duplicate SNPs are also removed here. Invalid non-crucial variables
  #	are set to missing.
  # The QC (in phases 2a, b and c) carries out three tests per variable:
  #	1) Does the column of this variable contain the correct (i.e.
  #	numeric or character) datatype? If not, the dataset is presumed
  #	corrupt and EmergencyExit is set to true.
  #	2) Which entries are missing?
  #	3) Which entries are invalid? (i.e. impossible values, like
  #		a p-value of 1.1 or chromosome 36)
  # The details of are different between variables, but this is the
  #	broad outline. This can make for a somewhat confusing nomenclature:
  #	-iNA-variables refer to data that is missing in the *original*
  #	dataset; while inv(alid) indicates entries with impossible
  #	values. inv-entries will later be set to NA (for effect allele, a 
  #	bad_ variable collects both invalid & missing entries).
  #	-_poor variables count the number of NA entries that have poor
  #	imputation quality and need not be reported as missing.
  
  
  # PHASE 2a: removing monomorphic SNPs
  #	& Y & M chromosome SNPs
  
  print("", quote = FALSE)
  print("Step 2: testing values and excluding SNPs", quote = FALSE)
  print(" - Checking alleles", quote = FALSE)
  flush.console()
  
  iNA_al1_list <- is.na(dataI$EFFECT_ALL)
  iNA_al1_L0	<- sum(iNA_al1_list)
  iNA_al2_list<-is.na(dataI$OTHER_ALL) | dataI$OTHER_ALL == "0" | dataI$OTHER_ALL == "-" | dataI$OTHER_ALL == "99"
  iNA_al2_N	<- sum(iNA_al2_list)
  
  if(iNA_al1_L0 == SNPn_input | iNA_al2_N == SNPn_input) {
    if(iNA_al1_L0 == SNPn_input) { 
      print("CRITICAL ERROR: effect-allele column is empty",quote=FALSE)
      save_log(2L, "allele data", "no effect allele", iNA_al1_L0, SNPn_input, "QC aborted", "Effect-allele column is empty.", filename_dir) 
    }
    if(iNA_al2_N == SNPn_input) { 
      print("CRITICAL ERROR: non-effect-allele column is empty",quote=FALSE)
      save_log(2L, "allele data", "no non-effect allele", iNA_al2_N, SNPn_input, "QC aborted", "Non-effect-allele column is empty.", filename_dir) 
    }
    EmergencyExit <- TRUE
  } else {
    if(!is.character(dataI$EFFECT_ALL)) {
      print("CRITICAL ERROR: effect-allele column contains non-character entries",quote=FALSE)
      save_log(2L, "allele data", "effect allele", SNPn_input - iNA_al1_L0, SNPn_input, "QC aborted", paste("Effect allele is", mode(dataI$EFFECT_ALL)), filename_dir)
      EmergencyExit <- TRUE
    }
    if(!is.character(dataI$OTHER_ALL)) {
      print("CRITICAL ERROR: non-effect-allele column contains non-character entries",quote=FALSE)
      save_log(2L, "allele data", "non-effect allele", SNPn_input - iNA_al2_N, SNPn_input, "QC aborted", paste("Non-effect allele is", mode(dataI$OTHER_ALL)), filename_dir)
      EmergencyExit <- TRUE
    }
  }
  
  # 1st phase-2 EmergencyExit: invalid or empty EFFECT_ALL/2 columns
  if(EmergencyExit) { return(emergency_return(name_in = filename_input, name_out = filename, LI = logI, LN = logN)) }
  
  # Removing spaces from alleles
  if(iNA_al1_L0 > 0L) {
    if(any(nchar(dataI$EFFECT_ALL) > 1L & !iNA_al1_list)) {
      dataI$EFFECT_ALL <- gsub("(^ +)|( +$)", "", dataI$EFFECT_ALL) }
  } else {
    if(any(nchar(dataI$EFFECT_ALL) > 1L)) {
      dataI$EFFECT_ALL <- gsub("(^ +)|( +$)", "", dataI$EFFECT_ALL) } }
  if(iNA_al2_N > 0L) {
    if(any(nchar(dataI$OTHER_ALL) > 1L & !iNA_al2_list)) {
      dataI$OTHER_ALL <- gsub("(^ +)|( +$)", "", dataI$OTHER_ALL) }
  } else {
    if(any(nchar(dataI$OTHER_ALL) > 1L)) {
      dataI$OTHER_ALL <- gsub("(^ +)|( +$)", "", dataI$OTHER_ALL) } }
  lwr_al1 <- which(dataI$EFFECT_ALL %in% c("a", "t", "c", "g"))
  lwr_al2 <- which(dataI$OTHER_ALL %in% c("a", "t", "c", "g"))
  if(length(lwr_al1) > 0L) { dataI$EFFECT_ALL[lwr_al1] <- toupper(dataI$EFFECT_ALL[lwr_al1])	}
  if(length(lwr_al2) > 0L) { dataI$OTHER_ALL[lwr_al2] <- toupper(dataI$OTHER_ALL[lwr_al2])	}
  
  inv_al2_list<- !dataI$OTHER_ALL %in% c("A", "T", "C", "G") & !iNA_al2_list
  inv_al2_N	<- sum(inv_al2_list)
  
  low_FRQ_list<- (dataI$EFF_ALL_FREQ == 0 | dataI$EFF_ALL_FREQ == 1) & !(is.na(dataI$EFF_ALL_FREQ) | iNA_al2_list | inv_al2_list)
  low_FRQ_N	<- sum(low_FRQ_list)
  
  same_al	<- which(dataI$EFFECT_ALL == dataI$OTHER_ALL & !(iNA_al2_list | inv_al2_list | low_FRQ_list))
  same_al_N	<- length(same_al)
  
  remove_L0	<- unique(c(which(iNA_al2_list | inv_al2_list | low_FRQ_list), same_al))
  monomorp_N <- length(remove_L0)
  
  if(is.character(dataI$CHR)) { dataI$CHR <- gsub("(^ +)|( +$)", "", dataI$CHR) }
  dataI$CHR[dataI$CHR %in% c("X", "x")] <- 23L
  dataI$CHR[dataI$CHR %in% c("Y", "y")] <- 24L
  dataI$CHR[dataI$CHR %in% c("XY","xy", "Xy", "xY")]<- 25L
  dataI$CHR[dataI$CHR %in% c("M", "m", "MT", "mt", "Mt", "mT")] <- 26L
  na_rm_ch <- any(is.na(dataI$CHR))
  
  if(remove_X) {
    chr_X_N <- sum(dataI$CHR == 23, na.rm = na_rm_ch)
    if(chr_X_N > 0L) {
      save_log(2L, "allele data", "Chromosome X SNPs", chr_X_N, SNPn_input, "SNPs removed", "-", filename_dir)
      print(paste(" - - Chromosome X SNPs removed:", chr_X_N), quote = FALSE)
      remove_L0 <- unique(c(remove_L0, which(dataI$CHR == 23)))
  } } else { chr_X_N <- NA }
  if(remove_Y) {
    chr_Y_N <- sum(dataI$CHR == 24, na.rm = na_rm_ch)
    if(chr_Y_N > 0L) {
      save_log(2L, "allele data", "Chromosome Y SNPs", chr_Y_N, SNPn_input, "SNPs removed", "-", filename_dir)
      print(paste(" - - Chromosome Y SNPs removed:", chr_Y_N), quote = FALSE)
      remove_L0 <- unique(c(remove_L0, which(dataI$CHR == 24)))
  } } else { chr_Y_N <- NA }
  if(remove_XY) {
    chr_XY_N <- sum(dataI$CHR == 25, na.rm = na_rm_ch)
    if(chr_XY_N > 0L) {
      save_log(2L, "allele data", "Chromosome XY SNPs", chr_XY_N, SNPn_input, "SNPs removed", "-", filename_dir)
      print(paste(" - - Chromosome XY SNPs removed:", chr_XY_N), quote = FALSE)
      remove_L0 <- unique(c(remove_L0, which(dataI$CHR == 25)))
  } } else { chr_XY_N <- NA }
  if(remove_M) {
    chr_M_N <- sum(dataI$CHR == 26, na.rm = na_rm_ch)
    if(chr_M_N > 0L) {
      save_log(2L, "allele data", "Chromosome M SNPs", chr_M_N, SNPn_input, "SNPs removed", "-", filename_dir)
      print(paste(" - - Chromosome M SNPs removed:", chr_M_N), quote = FALSE)
      remove_L0 <- unique(c(remove_L0, which(dataI$CHR == 26)))
  } } else { chr_M_N <- NA }
  
  remove_L0_N	<- length(remove_L0)
  
  # 2st step-2 EmergencyExit: no valid SNPs
  if(remove_L0_N > 0L) {
    print(paste(" - - monomorphic SNPs removed:", monomorp_N), quote = FALSE)
    if(iNA_al2_N > 0L) {
      save_log(2L, "allele data", "non-effect allele", iNA_al2_N, SNPn_input, "SNPs removed", "Missing non-effect allele entries", filename_dir)
      if(remove_L0_N == SNPn_input) { print(paste(" - - warning:", iNA_al2_N, "missing non-effect allele values"), quote = FALSE) }
    }
    if(inv_al2_N > 0L) {
      save_log(2L, "allele data", "non-effect allele", inv_al2_N, SNPn_input, "SNPs removed", "Invalid non-effect allele entries", filename_dir)
      print(paste(" - - warning:", inv_al2_N, "invalid non-effect allele values"), quote = FALSE)
      if(inv_al2_N > 30L) {	write.table(dataI[inv_al2_list, ][1:30, ]	, paste(filename_dir, "SNPs_invalid_OTHER_ALL.txt", sep = "_"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
      } else {			write.table(dataI[inv_al2_list, ]		, paste(filename_dir, "SNPs_invalid_OTHER_ALL.txt", sep = "_"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t") }
    }
    if(low_FRQ_N > 0L) {
      save_log(2L, "allele data", "allele frequency = 0", low_FRQ_N, SNPn_input, "SNPs removed", "Allele frequency of 0 or 1", filename_dir)
      if(remove_L0_N == SNPn_input) { print(paste(" - - warning:", low_FRQ_N, "SNPs with allele-frequency of 0 or 1"), quote = FALSE) }
    }
    if(same_al_N > 0L) {
      save_log(2L, "allele data", "identical alleles", same_al_N, SNPn_input, "SNPs removed", "Identical alleles, but allele frequency suggests non-monomorphic SNPs", filename_dir)
      print(paste(" - - warning:", same_al_N, "SNPs with identical alleles"), quote = FALSE)
    }
    if(remove_L0_N == SNPn_input) {
      print("CRITICAL ERROR: no usable SNPs in remaining dataset",quote=FALSE)
      save_log(2L, "allele data", "no usable SNPs", remove_L0_N, SNPn_input, "QC aborted", "No usable SNPs.", filename_dir)
      return(emergency_return(name_in = filename_input, name_out = filename, LI = logI, LN = logN))
    } else {
      save_log(2L, "allele data", "total removed", remove_L0_N, SNPn_input, "SNPs removed", "-", filename_dir)
      dataI <- dataI[-remove_L0, ]
    }
  }
  
  rm(iNA_al1_list, iNA_al2_list, lwr_al1, lwr_al2, inv_al2_list, low_FRQ_list, same_al, remove_L0)
  #	gc(verbose = FALSE) # memory cleaning
  

  #### Step 2b: checking imputation quality
  # inv and iNA are shadow-tables for dataI that log
  #	which entries are invalid / missing. Note that
  #	missing imputation-STATUS is tracked by inv, not iNA!
  # remove_L1_list will be filled in phase 2d. In 2b-c it
  # 	serves as template for the creation of other vectors.
  
  print(" - Checking data integrity", quote = FALSE)
  flush.console()
  
  SNPn_preQC		<- nrow(dataI)
  reason_excl		<- character(length = SNPn_preQC)
  column_improb	<- character(length = SNPn_preQC)
  remove_L1_list	<- logical(length = SNPn_preQC)
  
  inv <- data.frame(
    chr = remove_L1_list, pos = remove_L1_list,
    strand = remove_L1_list, p = remove_L1_list,
    FRQ = remove_L1_list, HWE = remove_L1_list,
    cal = remove_L1_list, N = remove_L1_list,
    impQ = remove_L1_list, impstatus = remove_L1_list )
  iNA <- data.frame(
    chr = is.na(dataI$CHR), pos = is.na(dataI$POSITION),
    strand = is.na(dataI$STRAND), p = is.na(dataI$PVALUE),
    FRQ = is.na(dataI$EFF_ALL_FREQ), HWE = is.na(dataI$HWE_PVAL),
    cal = is.na(dataI$CALLRATE), N = is.na(dataI$N_TOTAL),
    impQ = is.na(dataI$IMP_QUALITY))
  
  
  #Testing dataI$IMPUTED - unlike other parametes, the imputation status
  #	column can contain many data-types; hence it's "translated" here.
  #	The translated values will replace the original data at the
  #	end of phase 2. Imputation status is also unusual in that missing
  #	values are stored in the "invalid" table (because without
  #	this parameter call rate, HWE-p & imputation quality are useless), although
  #	the inv_impstatus_N still counts only true invalids. 
  if(is.character(dataI$IMPUTED)) {
    dataI$IMPUTED <- gsub("(^ +)|( +$)", "", dataI$IMPUTED) }
  iNA_impstatus_N <- if(length(imputed_NA) == 0L) 0L else sum(dataI$IMPUTED %in% imputed_NA)
  if(iNA_impstatus_N == SNPn_preQC) {
    save_log(2L, "data integrity", "no imputation status", iNA_impstatus_N, SNPn_preQC, "-", "Imputation-status column contains no values!", filename_dir)
    print(" - - warning: no imputation-status data", quote=FALSE)
    inv$impstatus <- TRUE
    column_improb <- paste0(column_improb, "imputation status; ")
    inv_impstatus_N <- 0L
  } else {
    new_impstatus <- convert_impstatus(dataI$IMPUTED, T_strings = imputed_T, F_strings = imputed_F, NA_strings = imputed_NA,
                                       use_log = TRUE, allSNPs = SNPn_preQC, fileL = filename_dir)
    inv$impstatus <- is.na(new_impstatus)
    column_improb[inv$impstatus] <- paste0(column_improb[inv$impstatus], "imputation status; ")
    inv_impstatus_N <- sum(inv$impstatus) - iNA_impstatus_N
    if(iNA_impstatus_N > 0L) {
      save_log(2L, "data integrity", "imputation status", iNA_impstatus_N, SNPn_preQC, "-", "Missing imputation status values", filename_dir)
      print(paste(" - - warning:", iNA_impstatus_N, "missing imputation-status values"), quote = FALSE)
    }
    if(inv_impstatus_N > 0L) {
      if(is.character(dataI$IMPUTED)) {
        print("CRITICAL ERROR: imputation-status column contains unidentified character strings",quote=FALSE)
        save_log(2L, "data integrity", "imputation status", inv_impstatus_N, SNPn_preQC, "QC aborted", "Unidentified character strings in imputation status column.", filename_dir)
        EmergencyExit <- TRUE
      } else {
        print(paste(" - - warning:", inv_impstatus_N, "invalid imputation-status values"), quote = FALSE)
        save_log(2L, "data integrity", "imputation status", inv_impstatus_N, SNPn_preQC, "set to NA", "Invalid imputation status", filename_dir)
  } } }
  if(inv_impstatus_N + iNA_impstatus_N == SNPn_preQC) {
    geno_list		<- remove_L1_list
    imp_list		<- remove_L1_list
    SNPn_preQC_geno	<- 0L
    SNPn_preQC_imp	<- 0L
  } else {
    geno_list		<- new_impstatus == 0 & !inv$impstatus
    imp_list		<- new_impstatus == 1 & !inv$impstatus
    SNPn_preQC_geno	<- sum(geno_list)
    SNPn_preQC_imp	<- sum( imp_list)
  }
  
  #Testing dataI$IMP_QUALITY
  iNA_impQ_N	 <- sum(iNA$impQ)
  if(!is.numeric(dataI$IMP_QUALITY) & iNA_impQ_N != SNPn_preQC) {
    EmergencyExit <- TRUE
    print("CRITICAL ERROR: imputation-quality column contains non-numeric entries",quote=FALSE)
    save_log(2L, "data integrity", "imputation quality", SNPn_preQC - iNA_impQ_N, SNPn_preQC, "QC aborted",	paste0("Imputation quality is ", mode(dataI$IMP_QUALITY)), filename_dir)
    low_impQ_list <- remove_L1_list
  } else {
    if(any(dataI$IMP_QUALITY == -1, na.rm = iNA_impQ_N > 0L)) {
      iNA$impQ[which(dataI$IMP_QUALITY == -1)] <- TRUE
      iNA_impQ_N <- sum(iNA$impQ)
      dataI$IMP_QUALITY[which(dataI$IMP_QUALITY == -1)] <- NA
    }
    iNA_impQ_g	<- sum(iNA$impQ & geno_list)
    iNA_impQ_i	<- sum(iNA$impQ &	imp_list)
    
    inv$impQ	<- !( (dataI$IMP_QUALITY >= minimal_impQ_value & dataI$IMP_QUALITY <= maximal_impQ_value) | iNA$impQ )
    inv_impQ_N  <- sum(inv$impQ)
    if(inv_impQ_N > 0L) { 
      column_improb[inv$impQ] <- paste0(column_improb[inv$impQ],"imputation quality; ")
      inv_impQ_g <- sum(inv$impQ & geno_list)
      inv_impQ_i <- sum(inv$impQ &	imp_list)
      save_log(2, "data integrity", "imputation quality", inv_impQ_N, SNPn_preQC, "set to NA", "Invalid imputation quality", filename_dir)
      print(paste(" - - warning:", inv_impQ_N, "invalid imputation-quality values"), quote = FALSE)
    } else {
      inv_impQ_g <- 0L
      inv_impQ_i <- 0L
    }
    if(ignore_impstatus) { low_impQ_list <- dataI$IMP_QUALITY < 0.3 & !iNA$impQ & !inv$impQ
    } else {               low_impQ_list <- dataI$IMP_QUALITY < 0.3 & !iNA$impQ & !inv$impQ & imp_list }
  }
  
  
  #### Step 2c: checking other variables
  
  # Testing dataI$MARKER for NA & duplicate entries.
  iNA_marker	<- which(is.na(dataI$MARKER))
  reason_excl[iNA_marker] <- paste(reason_excl[iNA_marker],"Missing marker ID;")
  iNA_mar_N	<- length(iNA_marker)
  if(!is.character(dataI$MARKER) & iNA_mar_N != SNPn_preQC) {
    EmergencyExit <- TRUE
    print("CRITICAL ERROR: marker-name column contains non-character entries",quote=FALSE)
    save_log(2L, "data integrity", "SNP names", SNPn_preQC, SNPn_preQC, "QC aborted", paste0("Marker name is ", mode(dataI$MARKER)), filename_dir)
  } else {
    dataI$MARKER <- gsub("(^ +)|( +$)", "", dataI$MARKER)
    dupli_list <- duplicated(dataI$MARKER)
    if(iNA_mar_N > 0L) {
      save_log(2L, "data integrity", "marker ID", iNA_mar_N, SNPn_preQC, "Markers removed", "Missing marker ID", filename_dir)
      print(paste(" - - warning:", iNA_mar_N, "missing marker names"), quote = FALSE)
      dupli_list[iNA_marker] <- FALSE	# removes "duplicate" NA's from dupli list
    }
    if(any(dupli_list)) {
      dupli		<- which(dataI$MARKER %in% unique(dataI$MARKER[dupli_list]) ) # Lists all SNPs with multiple entries
      reason_excl[dupli] <- paste(reason_excl[dupli],"Duplicate marker ID;")
      SNPn_preQC_dupli	<- length(dupli)
      save_log(2L, "data integrity", "duplicate marker IDs", SNPn_preQC_dupli, SNPn_preQC, "Markers removed", "Duplicate marker IDs", filename_dir)
      print(paste(" - - warning:", SNPn_preQC_dupli, "duplicate marker names"), quote = FALSE)
      write.table(dataI[dupli, ], paste(filename_dir, "SNPs_duplicates.txt", sep = "_"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
    } else {
      dupli	<- NULL
      SNPn_preQC_dupli <- 0L
    }
  }
  
  # Testing dataI$CHR
  iNA_chr_N <- sum(iNA$chr)
  if(iNA_chr_N != SNPn_preQC) {
    inv$chr	<- !(dataI$CHR %in% 1:26 | iNA$chr)
    column_improb[inv$chr] <- paste0(column_improb[inv$chr],"chromosome; ")
    inv_chr_N	<- sum(inv$chr)
    if(inv_chr_N > 0L) { 
      save_log(2L, "data integrity", "chromosome number", inv_chr_N, SNPn_preQC, "set to NA", "Invalid chromosome number", filename_dir)
      print(paste(" - - warning:", inv_chr_N, "invalid chromosome values"), quote = FALSE)
    }
    if(iNA_chr_N + inv_chr_N < SNPn_preQC) {
      if(any(!(c(1:22) %in% dataI$CHR))) {
        nochr <- paste(which(!(c(1:22) %in% dataI$CHR)), collapse=", ")
        print(paste(" - - no SNPs present for chromosome(s)", nochr), quote = FALSE)
        save_log(2L, "data integrity", "chromosome number", SNPn_preQC, SNPn_preQC, "-", paste("No SNPs present for chromosome(s)", nochr), filename_dir) 
    } }
  } else { inv_chr_N	<- 0L }
  
  #Testing dataI$POSITION
  iNA_pos_N	<- sum(iNA$pos)
  if(!is.numeric(dataI$POSITION) & iNA_pos_N != SNPn_preQC) {
    EmergencyExit <- TRUE
    print("CRITICAL ERROR: position column contains non-integer entries",quote=FALSE)
    save_log(2L, "data integrity", "position", SNPn_preQC, SNPn_preQC, "QC aborted", paste0("Position is ", mode(dataI$POSITION)), filename_dir)
  } else {
    inv$pos	 <- dataI$POSITION <= 0 & !iNA$pos
    column_improb[inv$pos] <- paste0(column_improb[inv$pos],"position; ")
    inv_pos_N <- sum(inv$pos)
    if(inv_pos_N > 0L) {
      save_log(2L, "data integrity", "position", inv_pos_N, SNPn_preQC, "set to NA", "Invalid positions", filename_dir)
      print(paste(" - - warning:", inv_pos_N, "invalid position values"), quote = FALSE)
  } }
  
  #Testing dataI$EFFECT_ALL
  # Because it is possible, if unlikely, that all non-missing entries were removed
  # in phase 2a, the data type is checked again, here.
  iNA_al1_N <- sum(is.na(dataI$EFFECT_ALL))
  if(!is.character(dataI$EFFECT_ALL) & iNA_al1_N != SNPn_preQC) {
    EmergencyExit <- TRUE
    print("CRITICAL ERROR: effect-allele column contains non-character entries",quote=FALSE)
    save_log(2L, "data integrity", "effect allele", SNPn_preQC, SNPn_preQC, "QC aborted", paste0("Effect allele is ", mode(dataI$EFFECT_ALL)), filename_dir)
  } else {
    bad_al1	<- which(!dataI[ ,"EFFECT_ALL"] %in% c("A", "T", "C", "G"))
    reason_excl[bad_al1] <- paste(reason_excl[bad_al1],"Invalid effect allele;")
    inv_al1_N	<- length(bad_al1) - iNA_al1_N
    if(iNA_al1_N > 0L) {
      save_log(2L, "data integrity", "effect allele", iNA_al1_N, SNPn_preQC, "Markers removed", "Missing effect allele", filename_dir)
      print(paste(" - - warning:", iNA_al1_N, "missing effect alleles"), quote = FALSE)
    }
    if(inv_al1_N > 0L) {
      save_log(2L, "data integrity", "effect allele", inv_al1_N, SNPn_preQC, "Markers removed", "Invalid effect allele", filename_dir)
      print(paste(" - - warning:", inv_al1_N, "invalid effect alleles"), quote = FALSE)      
    }
    
    if(!"A" %in% dataI$EFFECT_ALL) save_log(2L, "data integrity", "effect allele", SNPn_preQC, SNPn_preQC, "-", "No A alleles found in effect allele column", filename_dir)
    if(!"C" %in% dataI$EFFECT_ALL) save_log(2L, "data integrity", "effect allele", SNPn_preQC, SNPn_preQC, "-", "No C alleles found in effect allele column", filename_dir)
    if(!"T" %in% dataI$EFFECT_ALL) save_log(2L, "data integrity", "effect allele", SNPn_preQC, SNPn_preQC, "-", "No T alleles found in effect allele column", filename_dir)
    if(!"G" %in% dataI$EFFECT_ALL) save_log(2L, "data integrity", "effect allele", SNPn_preQC, SNPn_preQC, "-", "No G alleles found in effect allele column", filename_dir)
    
    if(!"A" %in% dataI$OTHER_ALL ) save_log(2L, "data integrity",  "other allele", SNPn_preQC, SNPn_preQC, "-", "No A alleles found in non-effect allele column", filename_dir)
    if(!"C" %in% dataI$OTHER_ALL ) save_log(2L, "data integrity",  "other allele", SNPn_preQC, SNPn_preQC, "-", "No C alleles found in non-effect allele column", filename_dir)
    if(!"T" %in% dataI$OTHER_ALL ) save_log(2L, "data integrity",  "other allele", SNPn_preQC, SNPn_preQC, "-", "No T alleles found in non-effect allele column", filename_dir)
    if(!"G" %in% dataI$OTHER_ALL ) save_log(2L, "data integrity",  "other allele", SNPn_preQC, SNPn_preQC, "-", "No G alleles found in non-effect allele column", filename_dir)
  }
  
  # Testing strand_info
  strand_pre <- list(plus = 0L, minus = 0L, missing = sum(iNA$strand), invalid = 0L)
  if(is.character(dataI$STRAND) & strand_pre$missing != SNPn_preQC) {
    dataI$STRAND <- gsub("(^ +)|( +$)", "", dataI$STRAND)
    if(strand_pre$missing > 0L) {
      strand_min <- dataI$STRAND == "-" & !iNA$strand
      inv$strand <- !(dataI$STRAND == "+" | strand_min | iNA$strand )
      column_improb[inv$strand] <- paste0(column_improb[inv$strand],"strand; ")
    } else {
      strand_min <- dataI$STRAND == "-"
      inv$strand <- !(dataI$STRAND == "+" | strand_min)
      column_improb[inv$strand] <- paste0(column_improb[inv$strand],"strand; ")
    }
    strand_pre$minus <- sum(strand_min)
    strand_pre$invalid <- sum(inv$strand)
    strand_pre$plus <- SNPn_preQC - (strand_pre$minus + strand_pre$missing + strand_pre$invalid)
    if(strand_pre$missing > 0L) { save_log(2L, "data integrity", "strand", strand_pre$missing, SNPn_preQC, "-", "Missing strand information", filename_dir) }
    if(strand_pre$invalid > 0L) {
      save_log(2L, "data integrity", "strand", strand_pre$invalid, SNPn_preQC, "set to NA", "Invalid strand information", filename_dir)
      print(paste(" - - warning:", strand_pre$invalid, "invalid strand values"), quote = FALSE)
    }
  } else {
    if(strand_pre$missing == SNPn_preQC) {
      save_log(2L, "data integrity", "strand", strand_pre$missing, SNPn_preQC, "-", "No strand information in dataset", filename_dir)
      print(" - - warning: no strand information in dataset", quote = FALSE)
      strand_min <- remove_L1_list
    } else {
      strand_pre$invalid <- SNPn_preQC - strand_pre$missing
      print("CRITICAL ERROR: strand column contains non-character entries", quote=FALSE)
      save_log(2L, "data integrity", "strand", strand_pre$invalid, SNPn_preQC, "QC aborted", paste0("Strand information is ", mode(dataI$STRAND)), filename_dir)
      EmergencyExit <- TRUE
  } }
  
  #Testing dataI$EFFECT
  iNA_eff	 <- which(is.na(dataI$EFFECT))
  reason_excl[iNA_eff] <- paste(reason_excl[iNA_eff],"Missing effect size;")
  iNA_eff_N <- length(iNA_eff)
  if(!is.numeric(dataI$EFFECT) & iNA_eff_N != SNPn_preQC) {
    EmergencyExit <- TRUE
    print("CRITICAL ERROR: effect size column contains non-numeric entries",quote=FALSE)
    save_log(2L, "data integrity", "effect size", SNPn_preQC, SNPn_preQC, "QC aborted", paste0("Effect size is ", mode(dataI$EFFECT)), filename_dir)
  } else {
    min_eff <- which(dataI$EFFECT == -1 & (dataI$PVALUE == -1 | dataI$STDERR == -1))
    if(length(min_eff) > 0L) {
      dataI$EFFECT[min_eff] <- NA
      iNA_eff	<- c(iNA_eff, min_eff)
      iNA_eff_N	<- iNA_eff_N + length(min_eff)
    }
    if(iNA_eff_N > 0L) {
      iNA_eff_poor	<- sum(low_impQ_list[iNA_eff])
      if(iNA_eff_N - iNA_eff_poor > 0L) {
        save_log(2L, "data integrity", "effect size", iNA_eff_N - iNA_eff_poor, SNPn_preQC, "Markers removed", "Missing effect size (poorly imputed markers not included)", filename_dir)
        print(paste(" - - warning:", iNA_eff_N - iNA_eff_poor, "missing effect size values"), quote = FALSE)
      }
    } else { iNA_eff_poor <- 0L }
  }
  
  #Testing dataI$STDERR
  iNA_se		<- which(is.na(dataI$STDERR))
  reason_excl[iNA_se] <- paste(reason_excl[iNA_se],"Missing standard error;")
  iNA_se_N <- length(iNA_se)
  if(!is.numeric(dataI$STDERR) & iNA_se_N != SNPn_preQC) {
    EmergencyExit <- TRUE
    print("CRITICAL ERROR: SE column contains non-numeric entries",quote=FALSE)
    save_log(2L, "data integrity", "standard error", SNPn_preQC, SNPn_preQC, "QC aborted", paste0("Standard error is ", mode(dataI$STDERR)), filename_dir)
  } else {
    if(any(dataI$STDERR == -1, na.rm = iNA_se_N > 0L)) {
      iNA_se						<- c(iNA_se, which(dataI$STDERR == -1))
      iNA_se_N						<- length(iNA_se)
      dataI$STDERR[which(dataI$STDERR == -1)]	<- NA
    }
    inv_se	<- which(dataI$STDERR <= 0)
    reason_excl[inv_se] <- paste(reason_excl[inv_se],"Invalid standard error;")
    inv_se_N	<- length(inv_se)			
    if(inv_se_N > 0L) {
      if(any(dataI$STDERR[inv_se] == 0)) {
        if(any(dataI$STDERR[inv_se] != 0)) {
          save_log(2L, "data integrity", "standard error", sum(dataI$STDERR[inv_se] != 0), SNPn_preQC, "Markers removed", "Invalid standard error (not counting SE = 0)", filename_dir)
        }
        save_log(2L, "data integrity", "standard error", sum(dataI$STDERR[inv_se] == 0), SNPn_preQC, "Markers removed", "Invalid standard error (SE = 0)", filename_dir)
      } else { 
        save_log(2L, "data integrity", "standard error", inv_se_N, SNPn_preQC, "Markers removed", "Invalid standard error", filename_dir)
      }
      print(paste(" - - warning:", inv_se_N, "invalid standard error values"), quote = FALSE)
    }
    if(iNA_se_N > 0L) {
      iNA_se_poor	<- sum(low_impQ_list[iNA_se])
      if(iNA_se_N - iNA_se_poor > 0L) {
        save_log(2L, "data integrity", "standard error", iNA_se_N - iNA_se_poor, SNPn_preQC, "Markers removed", "Missing standard error (poorly imputed markers not included)", filename_dir)
        print(paste(" - - warning:", iNA_se_N - iNA_se_poor, "missing standard error values"), quote = FALSE)
      }
    } else { iNA_se_poor<- 0L }
  }
  
  #Testing dataI$PVALUE
  iNA_p_N		<- sum(iNA$p)
  if(!is.numeric(dataI$PVALUE) & iNA_p_N != SNPn_preQC) {
    EmergencyExit <- TRUE
    print("CRITICAL ERROR: p-value column contains non-numeric entries",quote=FALSE)
    save_log(2L, "data integrity", "p-value", SNPn_preQC, SNPn_preQC, "QC aborted", paste0("P-value is ", mode(dataI$PVALUE)), filename_dir)
  } else {
    if(any(dataI$PVALUE == -1, na.rm = iNA_p_N > 0L)) {
      iNA$p[which(dataI$PVALUE == -1)]	<- TRUE
      iNA_p_N						<- sum(iNA$p)
      dataI$PVALUE[which(dataI$PVALUE == -1)]	<- NA
    }
    inv$p		<- !( (dataI$PVALUE > 0 & dataI$PVALUE <= 1) | iNA$p )
    column_improb[inv$p] <- paste0(column_improb[inv$p],"p-value; ")
    inv_p_N	<- sum(inv$p)
    if(inv_p_N > 0L) {
      save_log(2L, "data integrity", "p-value", inv_p_N, SNPn_preQC, "set to NA", "Invalid p-value", filename_dir)
      print(paste(" - - warning:", inv_p_N, "invalid p-values"), quote = FALSE)
  } }
  
  #Testing dataI$EFF_ALL_FREQ
  iNA_FRQ_N	<- sum(iNA$FRQ)
  if(!is.numeric(dataI$EFF_ALL_FREQ) & iNA_FRQ_N != SNPn_preQC) {
    EmergencyExit <- TRUE
    print("CRITICAL ERROR: allele-frequency column contains non-numeric entries",quote=FALSE)
    save_log(2L, "data integrity", "allele frequency", SNPn_preQC, SNPn_preQC, "QC aborted", paste0("Allele frequency is ", mode(dataI$EFF_ALL_FREQ)), filename_dir)
  } else {
    if(any(dataI$EFF_ALL_FREQ == -1, na.rm = iNA_FRQ_N > 0L)) {
      iNA$FRQ[which(dataI$EFF_ALL_FREQ == -1)] 	<- TRUE
      iNA_FRQ_N						<- sum(iNA$FRQ)
      dataI$EFF_ALL_FREQ[which(dataI$EFF_ALL_FREQ == -1)]	<- NA
    }
    inv$FRQ	<- !( (dataI$EFF_ALL_FREQ >= 0 & dataI$EFF_ALL_FREQ <= 1) | iNA$FRQ )
    column_improb[inv$FRQ] <- paste0(column_improb[inv$FRQ],"allele frequency; ")
    inv_FRQ_N	<- sum(inv$FRQ)
    if(inv_FRQ_N > 0L) {
      save_log(2L, "data integrity", "allele frequency", inv_FRQ_N, SNPn_preQC, "set to NA", "Invalid allele frequencies", filename_dir)
      print(paste(" - - warning:", inv_FRQ_N, "invalid allele frequencies"), quote = FALSE)
  } }
   
  #Testing dataI$HWE_PVAL
  iNA_HWE_N	<- sum(iNA$HWE)
  if(!is.numeric(dataI$HWE_PVAL) & iNA_HWE_N != SNPn_preQC) {
    EmergencyExit <- TRUE
    print("CRITICAL ERROR: HWE p-value column contains non-numeric entries",quote=FALSE)
    save_log(2L, "data integrity", "HWE p-value", SNPn_preQC, SNPn_preQC, "QC aborted", paste0("HWE p-value is ", mode(dataI$HWE_PVAL)), filename_dir)
  } else {
    if(any(dataI$HWE_PVAL == -1, na.rm = iNA_HWE_N > 0L)) {
      iNA$HWE[which(dataI$HWE_PVAL == -1)]		<- TRUE
      iNA_HWE_N							<- sum(iNA$HWE)
      dataI$HWE_PVAL[which(dataI$HWE_PVAL == -1)]	<- NA
    }
    inv$HWE <- !( (dataI$HWE_PVAL > 0 & dataI$HWE_PVAL <= 1) | iNA$HWE )
    column_improb[inv$HWE] <- paste0(column_improb[inv$HWE],"HWE p-value; ")
    inv_HWE_N <- sum(inv$HWE)
    
    iNA_HWE_g <- sum(iNA$HWE & geno_list)
    iNA_HWE_i <- sum(iNA$HWE & imp_list)
    if(inv_HWE_N > 0L) {
      save_log(2L, "data integrity", "HWE p-value", inv_HWE_N, SNPn_preQC, "set to NA", "Invalid HWE p-value", filename_dir)
      print(paste(" - - warning:", inv_HWE_N, "invalid HWE p-values"), quote = FALSE)
      inv_HWE_g <- sum(inv$HWE & geno_list)
      inv_HWE_i <- sum(inv$HWE & imp_list)
    } else {
      inv_HWE_g <- 0L
      inv_HWE_i <- 0L
  } }
  
  #Testing dataI$CALLRATE
  iNA_cal_N	<- sum(iNA$cal)
  if(!is.numeric(dataI$CALLRATE) & iNA_cal_N != SNPn_preQC) {
    EmergencyExit <- TRUE
    print("CRITICAL ERROR: call-rate column contains non-numeric entries",quote=FALSE)
    save_log(2L, "data integrity", "call rate", SNPn_preQC, SNPn_preQC, "QC aborted", paste0("Call rate is ", mode(dataI$CALLRATE)), filename_dir)
  } else {
    if(any(dataI$CALLRATE == -1, na.rm = iNA_cal_N > 0L)) {
      iNA$cal[which(dataI$CALLRATE == -1)]		<- TRUE
      iNA_cal_N							<- sum(iNA$cal)
      dataI$CALLRATE[which(dataI$CALLRATE == -1)]	<- NA
    }
    inv$cal	<- !( (dataI$CALLRATE >= 0 & dataI$CALLRATE <= 1) | iNA$cal )
    column_improb[inv$cal] <- paste0(column_improb[inv$cal],"call rate; ")
    inv_cal_N	<- sum(inv$cal)
    
    iNA_cal_g	<- sum(iNA$cal & geno_list)
    iNA_cal_i	<- sum(iNA$cal &	imp_list)
    if(inv_cal_N > 0L) {
      save_log(2L, "data integrity", "call rate", inv_cal_N, SNPn_preQC, "set to NA", "Invalid call rate", filename_dir)
      print(paste(" - - warning:", inv_cal_N, "invalid callrate values"), quote = FALSE)
      inv_cal_g <- sum(inv$cal & geno_list)
      inv_cal_i <- sum(inv$cal & imp_list)
    } else {
      inv_cal_g <- 0L
      inv_cal_i <- 0L
  } }

  #Testing dataI$N_TOTAL
  iNA_N_N	<- sum(iNA$N)
  if(!is.numeric(dataI$N_TOTAL) & iNA_N_N != SNPn_preQC) {
    EmergencyExit <- TRUE
    print("CRITICAL ERROR: QC aborted because sample-size column contains non-numeric entries",quote=FALSE)
    save_log(2L, "data integrity", "sample size", SNPn_preQC, SNPn_preQC, "QC aborted", paste0("Sample size is ", mode(dataI$N_TOTAL)), filename_dir)
  } else {
    iNA_N_g <- sum( iNA$N & geno_list)
    iNA_N_i <- sum( iNA$N & imp_list)
    inv$N	  <- dataI$N_TOTAL <= 0 & !iNA$N
    column_improb[inv$N] <- paste0(column_improb[inv$N],"sample size; ")
    inv_N_N <- sum(inv$N)
    if(inv_N_N > 0L) {
      save_log(2L, "data integrity", "sample size", inv_N_N, SNPn_preQC, "set to NA", "Invalid sample size", filename_dir)
      print(paste(" - - warning:", inv_N_N, "invalid sample size values"), quote = FALSE)
      inv_N_g <- sum(inv$N & geno_list)
      inv_N_i <- sum(inv$N &	imp_list)
    } else {
      inv_N_g <- 0L
      inv_N_i <- 0L
  } }
  
  # 3rd phase-2 EmergencyExit: invalid data-types
  if(EmergencyExit) return(emergency_return(name_in = filename_input, name_out = filename, LI = logI, LN = logN))
  
  
  #### PHASE 2d: removal of unusable SNPs
  #	Here the scripts checks if no column contained unexpected data
  #	(via EmergencyExit), and makes lists of entries selected for
  #	removal. These entries are then saved separately, together 
  #	with a list (inv_entries) of SNPs with improbable values for
  #	for non-crucial data. The improbable values are then set to NA,
  #	and the remove_L1 removed from the dataset.(In the
  #	unlikely event that all SNPs are removed, the script will
  #	activate EmergencyExit instead.)
  
  flush.console()
  SNPn_preQC_inv <- sum(rowSums(inv) > 0L)
  
  # 4rd phase-2 EmergencyExit: no usuable SNPs
  remove_L1 <- sort(unique(c(iNA_marker, dupli, bad_al1, iNA_eff, iNA_se, inv_se)))
  remove_L1_N <- length(remove_L1)
  if(remove_L1_N > 0L ) {
    if(remove_L1_N == SNPn_preQC) {
      print("CRITICAL ERROR: no usable markers remaining in dataset",quote=FALSE)
      save_log(2L, "data integrity", "unusable dataset", remove_L1_N, SNPn_preQC, "QC check aborted", "No valid markers.", filename_dir)
      return(emergency_return(name_in = filename_input, name_out = filename, LI = logI, LN = logN))
    } else {
      write.table(data.frame(dataI[remove_L1,],REASON=reason_excl[remove_L1]), paste(filename_dir, "SNPs_removed.txt", sep = "_"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
      
      if(SNPn_preQC_inv > 0L) {
        remove_L1_list[remove_L1] <- TRUE
        if(any(rowSums(inv) > 0L & !remove_L1_list)) {
          write.table(data.frame(dataI[rowSums(inv) > 0L & !remove_L1_list, ], COLUMN = column_improb[rowSums(inv) > 0L & !remove_L1_list], stringsAsFactors = FALSE)[1:min(1000, sum(rowSums(inv) > 0L & !remove_L1_list)), ],
                      paste(filename_dir, "SNPs_improbable_values.txt", sep = "_"), append = FALSE, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
        }
        
        if(inv_chr_N > 0L)          dataI$CHR[inv$chr] <- NA
        if(inv_pos_N > 0L)          dataI$POSITION[inv$pos] <- NA
        if(strand_pre$invalid > 0L) dataI$STRAND[inv$strand] <- NA
        if(inv_p_N > 0L)            dataI$PVALUE[inv$p] <- NA
        if(inv_FRQ_N > 0L)          dataI$EFF_ALL_FREQ[inv$FRQ] <- NA
        if(inv_HWE_N > 0L)          dataI$HWE_PVAL[inv$HWE] <- NA
        if(inv_cal_N > 0L)          dataI$CALLRATE[inv$cal] <- NA
        if(inv_N_N > 0L)            dataI$N_TOTAL[inv$N] <- NA
        if(inv_impQ_N > 0L)         dataI$IMP_QUALITY[inv$impQ] <- NA
      }
      dataI$IMPUTED <- if(iNA_impstatus_N == SNPn_preQC) NA else new_impstatus
            
      dataI	<- dataI[-remove_L1, ]
      inv <- inv[-remove_L1, ]
      iNA <- iNA[-remove_L1, ]
      geno_list <- geno_list[-remove_L1]
      imp_list <- imp_list[-remove_L1]
      strand_min <- strand_min[-remove_L1]
      SNPn_ref <- nrow(dataI)
      strand_ref <- list(
        plus = 0L,
        minus = if(strand_pre$minus == 0L) 0L else sum(strand_min),
        missing = if(strand_pre$missing == 0L) 0L else sum(iNA$strand),
        invalid = if(strand_pre$invalid == 0L) 0L else sum(inv$strand) )
      strand_ref$plus <- SNPn_ref - (strand_ref$missing + strand_ref$invalid + strand_ref$minus)
      save_log(2L, "data integrity", "total removed", remove_L1_N, SNPn_preQC, "Markers removed", "Total number of markers removed", filename_dir)
      print(paste(" - - unusable markers removed:", remove_L1_N), quote = FALSE)
    }
  } else {
    if(SNPn_preQC_inv > 0L) {
      write.table(data.frame(dataI[rowSums(inv) > 0L, ], COLUMN = column_improb[rowSums(inv) > 0L], stringsAsFactors = FALSE)[1:min(1000,SNPn_preQC_inv), ],
                  paste(filename_dir, "SNPs_improbable_values.txt", sep = "_"), append = FALSE, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
      
      if(inv_chr_N > 0L)          dataI$CHR[inv$chr] <- NA
      if(inv_pos_N > 0L)          dataI$POSITION[inv$pos] <- NA
      if(strand_pre$invalid > 0L) dataI$STRAND[inv$strand] <- NA
      if(inv_p_N > 0L)            dataI$PVALUE[inv$p] <- NA
      if(inv_FRQ_N > 0L)          dataI$EFF_ALL_FREQ[inv$FRQ] <- NA
      if(inv_HWE_N > 0L)          dataI$HWE_PVAL[inv$HWE] <- NA
      if(inv_cal_N > 0L)          dataI$CALLRATE[inv$cal] <- NA
      if(inv_N_N > 0L)            dataI$N_TOTAL[inv$N] <- NA
      if(inv_impQ_N > 0L)         dataI$IMP_QUALITY[inv$impQ] <- NA
    }
    dataI$IMPUTED <- if(iNA_impstatus_N == SNPn_preQC) NA else new_impstatus
    SNPn_ref <- SNPn_preQC
    strand_ref <- strand_pre
    save_log(2L, "data integrity", "total removed", 0L, SNPn_preQC, "-", "No markers removed", filename_dir)
    print(paste(" - - no markers removed from", filename_input), quote = FALSE)
  }
  
  if(iNA_impstatus_N != SNPn_preQC) rm(new_impstatus)
  rm(remove_L1, remove_L1_list, low_impQ_list,
     iNA_marker, dupli_list, dupli, bad_al1,
     iNA_eff, min_eff, iNA_se, inv_se, column_improb, reason_excl)
  gc(verbose = FALSE) # memory cleaning
  
  
  ### PHASE 3: checking alleles & frequency against a base (allele_ref_std) and
  #	a additional reference (allele_ref_alt) table. Base table will remain unaltered, while
  #	the additional table will be updated with unknown SNPs.
  
  print("", quote = FALSE)
  print("Step 3: checking alleles", quote = FALSE)
  
  # creating a high-quality list for the allele-frequency plots
  if(use_allele_std | use_allele_alt) {
    FRQ_FRQ_threshold <- if(useFRQ_threshold <= 1) SNPn_ref * useFRQ_threshold else useFRQ_threshold
    FRQ_HWE_threshold <- if(useHWE_threshold <= 1) SNPn_ref * useHWE_threshold else useHWE_threshold
    FRQ_cal_threshold <- if(useCal_threshold <= 1) SNPn_ref * useCal_threshold else useCal_threshold
    FRQ_imp_threshold <- if(useImp_threshold <= 1) SNPn_ref * useImp_threshold else useImp_threshold
    if(FRQ_FRQ_threshold < 2) {	FRQ_FRQ_threshold <- 2
                                save_log(3L, "thresholds", "allele frequency", SNPn_ref, SNPn_ref, "set to 2", "User-specified threshold was < 2", filename_dir) }
    if(FRQ_HWE_threshold < 2) {	FRQ_HWE_threshold <- 2
                                save_log(3L, "thresholds", "HWE p-value", SNPn_ref, SNPn_ref, "set to 2", "User-specified threshold was < 2", filename_dir) }
    if(FRQ_cal_threshold < 2) {	FRQ_cal_threshold <- 2
                                save_log(3L, "thresholds", "Call rate", SNPn_ref, SNPn_ref, "set to 2", "User-specified threshold was < 2", filename_dir) }
    if(FRQ_imp_threshold < 2) {	FRQ_imp_threshold <- 2
                                save_log(3L, "thresholds", "Imputation quality", SNPn_ref, SNPn_ref, "set to 2", "User-specified threshold was < 2", filename_dir) }
    if(ignore_impstatus) {    
      FRQ_HQ_list <- HQ_filter(data = dataI, ignore_impstatus = TRUE,
                               FRQ_val = if(SNPn_ref - sum(iNA$FRQ | inv$FRQ ) >= FRQ_FRQ_threshold) HQfilter_FRQ else NULL,
                               HWE_val = if(SNPn_ref - sum(iNA$HWE | inv$HWE ) >= FRQ_HWE_threshold) HQfilter_HWE else NULL,
                               cal_val = if(SNPn_ref - sum(iNA$cal | inv$cal ) >= FRQ_cal_threshold) HQfilter_cal else NULL,
                               imp_val = if(SNPn_ref - sum(iNA$impQ| inv$impQ) >= FRQ_imp_threshold) HQfilter_imp else NULL,
                               FRQ_NA = NAfilter_FRQ, HWE_NA = NAfilter_HWE, cal_NA = NAfilter_cal, imp_NA = NAfilter_imp)
    } else {
      FRQ_HQ_list <- HQ_filter(data = dataI, ignore_impstatus = FALSE,
                               FRQ_val = if(SNPn_ref - sum(iNA$FRQ | inv$FRQ) >= FRQ_FRQ_threshold) HQfilter_FRQ else NULL,
                               HWE_val = if(SNPn_preQC_geno == 0L) NULL else { if(sum(geno_list & !(iNA$HWE | inv$HWE) ) >= FRQ_HWE_threshold) HQfilter_HWE else NULL },
                               cal_val = if(SNPn_preQC_geno == 0L) NULL else { if(sum(geno_list & !(iNA$cal | inv$cal) ) >= FRQ_cal_threshold) HQfilter_cal else NULL },
                               imp_val = if(SNPn_preQC_imp  == 0L) NULL else { if(sum( imp_list & !(iNA$impQ|inv$impQ) ) >= FRQ_imp_threshold) HQfilter_imp else NULL },
                               FRQ_NA = NAfilter_FRQ, HWE_NA = NAfilter_HWE, cal_NA = NAfilter_cal, imp_NA = NAfilter_imp)
    }
  }
  
  # Step 3a: using the base reference
  if(use_allele_std) {
    ref_std <- dataI$MARKER %in% allele_ref_std$SNP
    SNPn_ref_std <- sum(ref_std)
    if(SNPn_ref_std > 0L) {
      print(paste(" - matching alleles to", allele_name_std), quote = FALSE)    	
      flush.console()
      allele_out_std <- match_alleles(dataset = dataI[ref_std,
                                                      c("MARKER", "EFFECT_ALL", "OTHER_ALL", "STRAND", "EFFECT", "EFF_ALL_FREQ")[c(TRUE, TRUE, TRUE, strand_ref$minus > 0L, TRUE, TRUE)]],
                                      dataname = filename,
                                      ref_set = allele_ref_std, ref_name = allele_name_std,
                                      unmatched_data = FALSE, HQ_subset = FRQ_HQ_list[ref_std],
                                      check_strand = strand_ref$minus > 0L,
                                      save_mismatches = TRUE, delete_mismatches = remove_mismatches_std,
                                      delete_diffEAF = remove_diffEAF_std, threshold_diffEAF = threshold_diffEAF,
                                      check_FRQ = TRUE, check_ambiguous = check_ambiguous_alleles,
                                      plot_FRQ = make_plots, plot_intensity = plot_intensity,
                                      plot_if_threshold = only_plot_if_threshold, threshold_r = threshold_allele_freq_correlation,
                                      return_SNPs = TRUE,
                                      save_name = filename, save_dir = dir_output, use_log = TRUE, log_SNPall = SNPn_ref )
      
      dataI$EFFECT_ALL[ref_std] <- allele_out_std$EFFECT_ALL
      dataI$OTHER_ALL[ref_std] <- allele_out_std$OTHER_ALL
      if(strand_ref$minus > 0L & allele_out_std$n_negative_strand > 0L) {
        dataI$STRAND[ref_std] <- allele_out_std$STRAND }
      dataI$EFFECT[ref_std] <- allele_out_std$EFFECT
      dataI$EFF_ALL_FREQ[ref_std] <- allele_out_std$EFF_ALL_FREQ
      
      # Memory-cleaning: data columns are removed from allele_out
      allele_out_std <- allele_out_std[1:16]
      gc(verbose = FALSE)
    } else {
      print(paste(" - matching alleles to", allele_name_std, "(SKIPPED)"), quote = FALSE)
      allele_out_std <- list(FRQ_cor = NA, FRQ_cor_ambiguous = NA, FRQ_cor_nonambi = NA,
                             n_SNPs = 0L, n_negative_strand = NA, n_negative_switch = NA, n_negative_mismatch = NA,
                             n_strandswitch = NA, n_mismatch = NA, n_flipped = NA, n_ambiguous = NA, n_suspect = NA, n_diffEAF = NA)
      save_log(3L, paste("allele match -", allele_name_std), "no data", SNPn_ref, SNPn_ref, "phase 3a skipped", paste("No matching markers found in", allele_name_std), filename_dir)
    }
  } else {
    SNPn_ref_std <- 0L
    allele_out_std <- list(FRQ_cor = NA, FRQ_cor_ambiguous = NA, FRQ_cor_nonambi = NA,
                           n_SNPs = NA, n_negative_strand = NA, n_negative_switch = NA, n_negative_mismatch = NA,
                           n_strandswitch = NA, n_mismatch = NA, n_flipped = NA, n_ambiguous = NA, n_suspect = NA, n_diffEAF = NA)
  }
  
  # Step 3b: using alternative reference
  if(use_allele_alt) {
    ref_alt <- if(use_allele_std) !ref_std & dataI$MARKER %in% allele_ref_alt$SNP else dataI$MARKER %in% allele_ref_alt$SNP
    SNPn_ref_alt<- sum(ref_alt)
    
    if(SNPn_ref_alt > 0L) {
      print(paste(" - matching alleles to", allele_name_alt), quote = FALSE)
      flush.console()
      allele_out_alt <- match_alleles(dataset = dataI[ref_alt, c("MARKER", "EFFECT_ALL", "OTHER_ALL", "STRAND", "EFFECT", "EFF_ALL_FREQ")[c(TRUE, TRUE, TRUE, strand_ref$minus > 0L, TRUE, TRUE)]],
                                      dataname = filename,
                                      ref_set = allele_ref_alt, ref_name = allele_name_alt,
                                      unmatched_data = FALSE, HQ_subset = FRQ_HQ_list[ref_alt],
                                      check_strand = strand_ref$minus > 0L,
                                      save_mismatches = TRUE, delete_mismatches = remove_mismatches_alt,
                                      delete_diffEAF = remove_diffEAF_alt, threshold_diffEAF = threshold_diffEAF,
                                      check_FRQ = TRUE, check_ambiguous = check_ambiguous_alleles,
                                      plot_FRQ = make_plots, plot_intensity = plot_intensity,
                                      plot_if_threshold = only_plot_if_threshold, threshold_r = threshold_allele_freq_correlation,
                                      return_SNPs = TRUE,
                                      save_name = filename, save_dir = dir_output, use_log = TRUE, log_SNPall = SNPn_ref )
      
      dataI$EFFECT_ALL[ref_alt] <- allele_out_alt$EFFECT_ALL
      dataI$OTHER_ALL[ref_alt] <- allele_out_alt$OTHER_ALL
      if(strand_ref$minus > 0L & allele_out_alt$n_negative_strand > 0L) { dataI$STRAND[ref_alt] <- allele_out_alt$STRAND }
      dataI$EFFECT[ref_alt] <- allele_out_alt$EFFECT
      dataI$EFF_ALL_FREQ[ref_alt]<- allele_out_alt$EFF_ALL_FREQ
      
      # Memory-cleaning: data columns are removed from allele_out
      allele_out_alt <- allele_out_alt[1:16]
      gc(verbose = FALSE)
    } else {
      print(paste(" - matching alleles to", allele_name_alt, "(SKIPPED)"), quote = FALSE)
      allele_out_alt <- list(FRQ_cor = NA, FRQ_cor_ambiguous = NA, FRQ_cor_nonambi = NA,
                             n_SNPs = 0L, n_negative_strand = NA, n_negative_switch = NA, n_negative_mismatch = NA,
                             n_strandswitch = NA, n_mismatch = NA, n_flipped = NA, n_ambiguous = NA, n_suspect = NA, n_diffEAF = NA)
    }
  } else {
    SNPn_ref_alt <- 0L
    allele_out_alt <- list(FRQ_cor = NA, FRQ_cor_ambiguous = NA, FRQ_cor_nonambi = NA,
                           n_SNPs = NA, n_negative_strand = NA, n_negative_switch = NA, n_negative_mismatch = NA,
                           n_strandswitch = NA, n_mismatch = NA, n_flipped = NA, n_ambiguous = NA, n_suspect = NA, n_diffEAF = NA)
  }

  # Step 3c: adding unknown SNPs to alternative reference
  SNPn_ref_new <- SNPn_ref - SNPn_ref_std - SNPn_ref_alt
  if(SNPn_ref_new > 0L) {
    if(update_alt) {
      if(use_allele_alt) {
        save_log(3L, paste("allele match -", allele_name_alt), "Markers not found", SNPn_ref_new, SNPn_ref, "added to alternative reference", paste("Markers not found in", allele_name_alt), filename_dir)
        print(paste(" - adding new SNPs to", allele_name_alt), quote = FALSE)
      } else {
        save_log(3L, paste("allele match -", allele_name_alt), "Markers not found", SNPn_ref_new, SNPn_ref, "created alternative reference", paste("Markers not found: creating", allele_name_alt), filename_dir)
        print(paste(" - creating", allele_name_alt, "file"), quote = FALSE)
      }
    } else if (use_allele_alt) {
      print(" - checking remaining alleles", quote = FALSE)
    }
    flush.console()
    
    if(SNPn_ref_std > 0L | SNPn_ref_alt > 0L) {
      if(SNPn_ref_std > 0L & SNPn_ref_alt > 0L) {
        ref_new <- !(ref_std | ref_alt)
      } else {
        if(SNPn_ref_std > 0L) {	ref_new <- !ref_std
        } else {			ref_new <- !ref_alt }
      }
    } else {					ref_new <- !logical(length = SNPn_ref) }
    
    allele_out_new <- list(n_SNPs = SNPn_ref_new,
                           n_negative_strand = if(strand_ref$minus == 0L) { 0L } else { sum(strand_min[ref_new]) },
                           n_flipped = sum(dataI$EFF_ALL_FREQ[ref_new] > 0.5, na.rm = iNA_FRQ_N+inv_FRQ_N>0L),
                           n_ambiguous = sum(( dataI$EFFECT_ALL[ref_new] == "A" & dataI$OTHER_ALL[ref_new] == "T" ) | ( dataI$EFFECT_ALL[ref_new] == "T" & dataI$OTHER_ALL[ref_new] == "A" ) | ( dataI$EFFECT_ALL[ref_new] == "C" & dataI$OTHER_ALL[ref_new] == "G" ) | ( dataI$EFFECT_ALL[ref_new] == "G" & dataI$OTHER_ALL[ref_new] == "C" ) ) )
    
    if(allele_out_new$n_negative_strand > 0L) {
      dataI[ref_new & strand_min, c("EFFECT_ALL", "OTHER_ALL", "STRAND")] <-
        switch_strand(dataI[ref_new & strand_min, c("EFFECT_ALL", "OTHER_ALL", "STRAND")], strand_col = TRUE)
    }
    
    if(allele_out_new$n_flipped > 0L) {
      ref_new_flip <- ref_new & dataI$EFF_ALL_FREQ > 0.5 & !iNA$FRQ & !inv$FRQ
      ref_new_al2	 <- dataI[ref_new_flip, "OTHER_ALL"]
      dataI[ref_new_flip, "OTHER_ALL"] <- dataI[ref_new_flip, "EFFECT_ALL"]
      dataI[ref_new_flip, "EFFECT_ALL"] <- ref_new_al2
      dataI[ref_new_flip, "EFFECT"]	<- -dataI[ref_new_flip, "EFFECT"]
      dataI[ref_new_flip, "EFF_ALL_FREQ"] <- 1 - dataI[ref_new_flip, "EFF_ALL_FREQ"]
    }
    
    if(update_alt) {
      if(update_as_rdata) {
        update_savename <- paste0(update_savename, ".RData")
        if(backup_alt & file.exists(paste(dir_references, update_savename, sep = "/"))) {
          save_log(3L, "allele match ", "back-up", 0L, 1L, "saved back-up", paste0("Back-up of previous alternative reference file saved as: ", dir_references, "/", remove_extension(update_savename), "_", Sys.Date(), ".RData"), filename_dir)
          print(" - - creating back-up of previous alternative-reference file", quote = FALSE)
          flush.console()
          file.rename(paste(dir_references, update_savename, sep = "/"), paste0(dir_references, "/", remove_extension(update_savename), "_", Sys.Date(), ".RData"))
        }
        print(" - - saving alternative-reference file", quote = FALSE)
        flush.console()
        if(use_allele_alt) {
          allele_ref_alt <- rbind(allele_ref_alt, data.frame(SNP = dataI$MARKER[ref_new], CHR = dataI$CHR[ref_new], POS = as.integer(dataI$POSITION[ref_new]),
                                                             MINOR = as.factor(dataI$EFFECT_ALL[ref_new]), MAJOR = as.factor(dataI$OTHER_ALL[ref_new]), MAF = dataI$EFF_ALL_FREQ[ref_new],
                                                             SOURCE = as.factor(filename_input), DATE_ADDED = as.factor(date()), stringsAsFactors = FALSE) )
        } else {
          allele_ref_alt <- data.frame(SNP = dataI$MARKER[ref_new], CHR = dataI$CHR[ref_new], POS = as.integer(dataI$POSITION[ref_new]),
                                       MINOR = as.factor(dataI$EFFECT_ALL[ref_new]), MAJOR = as.factor(dataI$OTHER_ALL[ref_new]), MAF = dataI$EFF_ALL_FREQ[ref_new],
                                       SOURCE = as.factor(filename_input), DATE_ADDED = as.factor(date()), stringsAsFactors = FALSE)
        }
        save(allele_ref_alt, file = paste(dir_references, update_savename, sep = "/"))
        rm(allele_ref_alt)
        # In order save the correctobject name in the .RData file, it is
        # necesarry to edit allele_ref_alt itself. Since allelel_ref_alt
        # is often a pointer to an object outside of QC_GWAS, the object
        # is duplicated in the memory (one unaltered instance outside
        # QC_GWAS and one updated instance inside). This may use considerable
        # amounts of memory; hence allele_ref_alt is deleted once it is no longer
        # necessary.
        
      } else {
        update_savename <- paste0(update_savename, ".txt")
        if(backup_alt & file.exists(paste(dir_references, update_savename, sep = "/"))) {
          save_log(3L, "allele match ", "back-up", 1L, 1L, "created back-up", paste0("copy of previous alternative reference file saved as: ", dir_references, "/", remove_extension(update_savename), "_", Sys.Date(), ".RData"), filename_dir)
          print(" - - creating back-up of previous alternative-reference file", quote = FALSE)
          flush.console()
          file.copy(paste(dir_references, update_savename, sep = "/"), paste0(dir_references, "/", remove_extension(update_savename), "_", Sys.Date(), ".txt"))
        }
        print(" - - saving alternative-reference file", quote = FALSE)
        flush.console()
        if(use_allele_alt) {
          write.table(rbind(allele_ref_alt, data.frame(SNP = dataI$MARKER[ref_new], CHR = dataI$CHR[ref_new], POS = dataI$POSITION[ref_new],
                                                       MINOR = dataI$EFFECT_ALL[ref_new], MAJOR = dataI$OTHER_ALL[ref_new], MAF = dataI$EFF_ALL_FREQ[ref_new],
                                                       SOURCE = filename_input, DATE_ADDED = date(), stringsAsFactors = FALSE) ),
                      paste(dir_references, update_savename, sep ="/"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
        } else {
          write.table(data.frame(SNP = dataI$MARKER[ref_new], CHR = dataI$CHR[ref_new], POS = dataI$POSITION[ref_new],
                                 MINOR = dataI$EFFECT_ALL[ref_new], MAJOR = dataI$OTHER_ALL[ref_new], MAF = dataI$EFF_ALL_FREQ[ref_new],
                                 SOURCE = filename_input, DATE_ADDED = date(), stringsAsFactors = FALSE),
                      paste(dir_references, update_savename, sep ="/"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
        }
      }
      allele_ref_changed	<- TRUE
    }
  } else { allele_out_new <- list(n_SNPs = 0L, n_negative_strand = NA, n_flipped = NA, n_ambiguous = NA ) }
  
  #Step 3d: removal of mismatching and different allele frequency entries from 3a & 3b
  if(use_allele_std | use_allele_alt) {
    N_mismatch <- sum(allele_out_std$n_mismatch, allele_out_alt$n_mismatch, na.rm = !(use_allele_std & use_allele_alt))
    N_diffEAF <- sum(allele_out_std$n_diffEAF, allele_out_alt$n_diffEAF, na.rm = !(use_allele_std & use_allele_alt))
    N_ambiguous <- sum(allele_out_std$n_ambiguous, allele_out_alt$n_ambiguous, allele_out_new$n_ambiguous, na.rm = TRUE)
    N_suspect <- sum(allele_out_std$n_suspect, allele_out_alt$n_suspect, na.rm = TRUE)		
    if(N_suspect > 0L) save_log(3L, "allele match", "suspicious strand", N_suspect, SNPn_ref, "-", "Possible strand mismatch: ambiguous markers with minor allele frequency > 35%", filename_dir)
    
    use_L2 <- FALSE
    if(SNPn_ref_std > 0L & allele_out_std$n_mismatch > 0L) {
      print(paste0(" - - mismatches with ", allele_name_std, ": ", allele_out_std$n_mismatch), quote = FALSE)
      if(remove_mismatches_std) use_L2 <- TRUE
    }
    if(SNPn_ref_alt > 0L & allele_out_alt$n_mismatch > 0L) {
      print(paste0(" - - mismatches with ", allele_name_alt, ": ", allele_out_alt$n_mismatch), quote = FALSE)
      if(remove_mismatches_alt) use_L2 <- TRUE
    }
    if(SNPn_ref_std > 0L & allele_out_std$n_diffEAF > 0L) {
      print(paste0(" - - differing allele frequencies from ", allele_name_std, ": ", allele_out_std$n_diffEAF), quote = FALSE)
      if(remove_diffEAF_std) use_L2 <- TRUE
    }
    if(SNPn_ref_alt > 0L & allele_out_alt$n_diffEAF > 0L) {
      print(paste0(" - - differing allele frequencies from ", allele_name_alt, ": ", allele_out_alt$n_diffEAF), quote = FALSE)
      if(remove_diffEAF_alt) use_L2 <- TRUE
    }
    flush.console()
    
    if(use_L2) {
      remove_L2_list <- is.na(dataI$EFFECT_ALL)
      remove_L2_N <- sum(remove_L2_list)
      
      if(remove_L2_N == SNPn_ref) {
        # Step 3 EmergencyExit: all SNPs removed
        print("CRITICAL ERROR: no SNPs remaining in dataset", quote = FALSE)
        save_log(3L, "allele match", "no matching SNPs", SNPn_ref, SNPn_ref, "QC check aborted", "All markers excluded due to mismatching alleles or allele-frequencies", filename_dir)
        return(emergency_return(name_in = filename_input, name_out = filename, LI = logI, LN = logN))
      }
      
      dataI <- dataI[!remove_L2_list, ]
      inv <- inv[!remove_L2_list, ]
      iNA <- iNA[!remove_L2_list, ]
      geno_list <- geno_list[!remove_L2_list]
      imp_list  <- imp_list[!remove_L2_list]
      strand_min <- strand_min[!remove_L2_list]
      rm(remove_L2_list)   
    } else { remove_L2_N <- 0L }
    
    rm(FRQ_HQ_list)
    if(use_allele_std) rm(ref_std)
    if(use_allele_alt) rm(ref_alt)
  } else {
    N_mismatch <- NA
    N_diffEAF <- NA
    N_ambiguous <- allele_out_new$n_ambiguous
    N_suspect <- NA
    remove_L2_N <- 0L
  }
  
  if(SNPn_ref_new > 0L) rm(ref_new)
  gc(verbose = FALSE) # memory cleaning
  
  #### Step 4: QC of the final dataset
  # Step 4a
  # Determining whether sufficient data is available for the following QC tests
  print("", quote = FALSE)
  print("Step 4: generating QC statistics & graphs", quote = FALSE)
  flush.console()
  
  SNPn_postQC	<- nrow(dataI)
  SNPn_postQC_inv	<- sum(rowSums(inv) > 0L)
  SNPn_postQC_geno<- sum(geno_list)
  SNPn_postQC_imp	<- sum( imp_list)
  
  if(remove_L2_N > 0L) {
    strand_post <- list(
      plus = 0L,
      minus = if(strand_ref$minus == 0L) 0L else sum(strand_min),
      missing = if(strand_ref$missing == 0L) 0L else sum(iNA$strand),
      invalid = if(strand_ref$invalid == 0L) 0L else sum(inv$strand) )
    strand_post$plus <- SNPn_postQC - (strand_post$missing + strand_post$invalid + strand_post$minus)
  } else { strand_post <- strand_ref }
  
  if(useMan_threshold <= 1) {	useMan_threshold <- SNPn_postQC * useMan_threshold }
  if(useMan_threshold < 10) {	useMan_threshold <- 10
                              save_log(4L, "thresholds", "Manhattan", SNPn_postQC, SNPn_postQC, "set to 10", "User-specified threshold was < 10", filename_dir) }
  # Threshold is set to 10, so it corresponds to the QC_plots internal threshold
  if(useFRQ_threshold <= 1) {	useFRQ_threshold <- SNPn_postQC * useFRQ_threshold }
  if(useFRQ_threshold <  2) {	useFRQ_threshold <- 2
                             save_log(4L, "thresholds", "allele frequency", SNPn_postQC, SNPn_postQC, "set to 2", "User-specified threshold was < 2", filename_dir) }
  if(useHWE_threshold <= 1) {	useHWE_threshold <- SNPn_postQC * useHWE_threshold }
  if(useHWE_threshold <  2) {	useHWE_threshold <- 2
                             save_log(4L, "thresholds", "HWE p-value", SNPn_postQC, SNPn_postQC, "set to 2", "User-specified threshold was < 2", filename_dir) }
  if(useCal_threshold <= 1) {	useCal_threshold <- SNPn_postQC * useCal_threshold }
  if(useCal_threshold <  2) {	useCal_threshold <- 2
                             save_log(4L, "thresholds", "Call rate", SNPn_postQC, SNPn_postQC, "set to 2", "User-specified threshold was < 2", filename_dir) }
  if(useImp_threshold <= 1) {	useImp_threshold <- SNPn_postQC * useImp_threshold }
  if(useImp_threshold <  2) {	useImp_threshold <- 2
                              save_log(4L, "thresholds", "Imputation quality", SNPn_postQC, SNPn_postQC, "set to 2", "User-specified threshold was < 2", filename_dir) }
  
  useMan <- SNPn_postQC - sum(iNA$chr | inv$chr | iNA$pos | inv$pos) >= useMan_threshold
  useFRQ <- SNPn_postQC - sum(iNA$FRQ | inv$FRQ) >= useFRQ_threshold
  if(ignore_impstatus) {
    useHWE <- SNPn_postQC - sum(iNA$HWE | inv$HWE ) >= useHWE_threshold
    useCal <- SNPn_postQC - sum(iNA$cal | inv$cal ) >= useCal_threshold
    useImp <- SNPn_postQC - sum(iNA$impQ| inv$impQ) >= useImp_threshold
  } else {
    if(SNPn_postQC_geno == 0L) {
      useHWE <- FALSE
      useCal <- FALSE
    } else {
      useHWE <- sum(geno_list & !(iNA$HWE | inv$HWE) ) >= useHWE_threshold
      useCal <- sum(geno_list & !(iNA$cal | inv$cal) ) >= useCal_threshold
    }
    if(SNPn_postQC_imp == 0L) {
      useImp <- FALSE
    } else {
      useImp <- sum( imp_list & !(iNA$impQ|inv$impQ) ) >= useImp_threshold
  } }
  
  if(!useFRQ) {
    QQfilter_FRQ <- NULL
    HQfilter_FRQ <- NULL }
  if(!useHWE) {
    QQfilter_HWE <- NULL
    HQfilter_HWE <- NULL }
  if(!useCal) {
    QQfilter_cal <- NULL
    HQfilter_cal <- NULL }
  if(!useImp) {
    QQfilter_imp <- NULL
    HQfilter_imp <- NULL }
  
  HQ_list <- HQ_filter(data = dataI, ignore_impstatus = ignore_impstatus,
                       FRQ_val = HQfilter_FRQ, HWE_val = HQfilter_HWE,
                       cal_val = HQfilter_cal, imp_val = HQfilter_imp,
                       FRQ_NA = NAfilter_FRQ, HWE_NA = NAfilter_HWE, cal_NA = NAfilter_cal, imp_NA = NAfilter_imp)
  SNPn_postQC_HQ <- sum(HQ_list)
  useHQ	<- SNPn_postQC_HQ > 9L
  
  if(!useHQ) {
    save_log(4L, "QC statistics", "insufficient high-quality markers", SNPn_postQC, SNPn_postQC, "-", "Less than 10 high-quality SNPs: cannot apply standard filter", filename_dir)
    print(" - - warning: insufficient high-quality SNPs to apply HQ filter", quote = FALSE)
  }
  
  # Step 4b: histograms
  if(plot_histograms) {
    print(" - creating histograms", quote = FALSE)
    flush.console()
    
    png(paste(filename_dir, "graph_histogram.png", sep = "_"),	width = 1440, height = 960, res = 108)
    par(mfrow = c(2, 3 ))
    
    hist(dataI$EFFECT,
         breaks = 30, freq = FALSE, col = "blue", plot = TRUE,
         main = "Effect size", xlab = "Effect size", cex.main = 1.5)
    
    hist(dataI$STDERR,
         breaks = 30, freq = FALSE, col = "yellow", plot = TRUE,
         main = "Standard error", xlab = "Standard error", cex.main = 1.5)
    
    if(useFRQ) {
      hist(dataI$EFF_ALL_FREQ, breaks = 20, freq = FALSE, col = "green3", xlim = c(0,1), plot = TRUE,
           main = "Allele frequency", xlab = "Allele frequency", cex.main = 1.5)
    } else {
      plot(0.5, 0.5, xlim = c(0,1), ylim = c(0,1), col = "white", cex.main = 1.5,
           main = "Allele frequency", xlab = "Allele frequency", ylab = "Density")
      text(0.5, 0.5, "Insufficient data", cex = 1.5)
    }
    
    if(useHWE) {
      if(ignore_impstatus) {
        hist(dataI$HWE_PVAL, breaks = 20, freq = FALSE, col = "yellow", xlim = c(0,1), plot = TRUE,
             main = "HWE p-value", xlab = "HWE p-value", cex.main = 1.5)
        mtext("Genotyped & imputed SNPs", cex = 0.8, line = 0.5, col = "red")
      } else {
        hist(dataI$HWE_PVAL[geno_list], breaks = 20, freq = FALSE, col = "yellow", xlim = c(0,1), plot = TRUE,
             main = "HWE-p distribution", xlab = "HWE p-value", cex.main = 1.5)
        mtext("Genotyped SNPs only", cex = 0.8, line = 0.5, col = "blue") }
    } else {
      plot(0.5, 0.5 , xlim = c(0,1), ylim = c(0,1), col = "white", cex.main = 1.5,
           main = "HWE-p distribution", xlab = "HWE p-value", ylab = "Density")
      text(0.5, 0.5, "Insufficient data", cex = 1.5) }
    
    if(useCal) {
      if(ignore_impstatus) {
        hist(dataI$CALLRATE, freq = FALSE, col = "red", xlim = c(0.75,1), plot = TRUE,
             breaks = c(0, 0.8, 0.85, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1),
             main = "Call rate", xlab = "Call rate", cex.main = 1.5)
        mtext("Genotyped & imputed SNPs", cex = 0.8, line = 0.5, col = "red")
        if(min(dataI$CALLRATE, na.rm = iNA_cal_N + inv_cal_N > 0) < 0.75) {
          title(sub = paste("Note:", sum(dataI$CALLRATE < 0.75, na.rm = iNA_cal_N + inv_cal_N > 0), "outliers excluded from graph"), col.sub = "blue", cex.sub = 1.3) }
      } else {
        hist(dataI$CALLRATE[geno_list], freq = FALSE, col = "red", xlim = c(0.75,1), plot = TRUE,
             breaks = c(0, 0.8, 0.85, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1),
             main = "Call rate", xlab = "Call rate", cex.main = 1.5)
        mtext("Genotyped SNPs only", cex = 0.8, line = 0.5, col = "blue")
        if(min(dataI$CALLRATE[geno_list], na.rm = iNA_cal_g + inv_cal_g > 0) < 0.75) {
          title(sub = paste("Note:", sum(dataI$CALLRATE[geno_list] < 0.75, na.rm = iNA_cal_g + inv_cal_g > 0), "outliers excluded from graph"), col.sub = "blue", cex.sub = 1.3) } }
    } else {
      plot(0.5, 0.5, xlim = c(0,1), ylim = c(0,1), col = "white", cex.main = 1.5,
           main = "Call rate", xlab = "Call rate", ylab = "Density")
      text(0.5, 0.5, "Insufficient data", cex = 1.5) }
    
    if(useImp) {
      if(ignore_impstatus) {
        hist(dataI$IMP_QUALITY, freq = FALSE, col = "blue", xlim = c(floor(min(0, dataI$IMP_QUALITY, na.rm = iNA_impQ_N+inv_impQ_N>0L) * 10)/10, ceiling(max(dataI$IMP_QUALITY, na.rm = iNA_impQ_N+inv_impQ_N>0L) * 10)/10), plot = TRUE,
             main = "Imputation Quality", xlab = "Imputation quality", cex.main = 1.5)
        mtext("Genotyped & imputed SNPs", cex = 0.8, line = 0.5, col = "red")
      } else {
        hist(dataI$IMP_QUALITY[imp_list], freq = FALSE, col = "blue", xlim = c(floor(min(0, dataI$IMP_QUALITY[imp_list], na.rm = iNA_impQ_i+inv_impQ_i>0L) * 10)/10, ceiling(max(dataI$IMP_QUALITY[imp_list], na.rm = iNA_impQ_i+inv_impQ_i>0L) * 10)/10), plot = TRUE,
             main = "Imputation Quality", xlab = "Imputation quality", cex.main = 1.5)
        mtext("Imputed SNPs only", cex = 0.8, line = 0.5, col = "blue") }
    } else {
      plot(0.5, 0.5, xlim = c(0,1), ylim = c(0,1), col = "white", cex.main = 1.5, 
           main = "Imputation Quality", xlab = "Imputation quality", ylab = "Density")
      text(0.5, 0.5, "Insufficient data", cex = 1.5) }
    
    title(sub = filename, cex.sub = 1.3)
    dev.off()
    gc(verbose = FALSE)
  }
  
  outcome_ES <- list(HQ_SNPs = SNPn_postQC_HQ, p = NA, return_ES = return_HQ_effectsizes,
                     HQ_effectsizes = if(return_HQ_effectsizes) {
                       if(SNPn_postQC_HQ < 1000L) { c(dataI$EFFECT[HQ_list], rep(NA, 1000L - SNPn_postQC_HQ))
                       } else { sort((dataI$EFFECT[HQ_list])[c(1:1000)*(SNPn_postQC_HQ %/% 1000)]) } } else { NULL } )
  
  
  # Step 4c: checking p-values
  print(" - checking p-values", quote = FALSE)
  outcome_P <- check_P(dataI, threshold_r = threshold_p_correlation, plot_correlation = make_plots, plot_if_threshold = only_plot_if_threshold,
                       save_name = filename, save_dir = dir_output, dataN = SNPn_postQC, use_log = TRUE, HQ_subset = HQ_list)
  
  if(calculate_missing_p) {
    calc_p_newN <- sum(iNA$p | inv$p)
    if(calc_p_newN > 0L) {
      print(" - calculating missing/invalid p-values", quote = FALSE)
      dataI$PVALUE <- ifelse(iNA$p | inv$p, pchisq((dataI$EFFECT/dataI$STDERR)^2, 1, lower.tail=FALSE), dataI$PVALUE)
      save_log(4L, "p-value correction", "missing/invalid p-values", calc_p_newN, SNPn_postQC, "p-values recalculated", "Missing/invalid p-values replaced with recalculated p-values.", filename_dir)
      
      low_p_newN <- sum(dataI$PVALUE < 1e-300 & (iNA$p | inv$p) )
      if(low_p_newN > 0L) {
        dataI$PVALUE[dataI$PVALUE < 1e-300 & (iNA$p | inv$p) ] <- 1e-300
        save_log(4L, "p-value correction", "extreme p-values", low_p_newN, SNPn_postQC, "p-value set to 1e-300", "Recalculation resulted in extreme p-value: it has been set to 1e-300", filename_dir)
      }
    } else { low_p_newN <- 0L }
  } else {
    calc_p_newN <- NA
    low_p_newN <- 0L
  }
  
  # Step 4d: creating QQ & Manhattan plots
  plot_output <- QC_plots(dataset = dataI, plot_QQ = plot_QQ, plot_Man = useMan & plot_Manhattan,
                          FRQfilter_values = QQfilter_FRQ, FRQfilter_NA = NAfilter_FRQ,
                          HWEfilter_values = QQfilter_HWE, HWEfilter_NA = NAfilter_HWE,
                          calfilter_values = QQfilter_cal, calfilter_NA = NAfilter_cal,
                          impfilter_values = QQfilter_imp, impfilter_NA = NAfilter_imp, impfilter_min = minimal_impQ_value,
                          manfilter_FRQ = HQfilter_FRQ, manfilter_HWE = HQfilter_HWE, manfilter_cal = HQfilter_cal, manfilter_imp = HQfilter_imp,
                          ignore_impstatus = ignore_impstatus,
                          plot_QQ_bands = plot_QQ_bands, plot_cutoff_p = plot_cutoff_p, plot_names = FALSE,
                          save_name = filename, save_dir = dir_output, use_log = TRUE)
  gc(verbose = FALSE) # memory cleaning
  if(plot_output$lambda[1] > 1.1 & !is.na(plot_output$lambda[1])) {
    save_log(phaseL = 5L, checkL = "QC statistics", typeL = "high lambda", SNPL = sum(!is.na(dataI$PVALUE)), allSNPs = SNPn_postQC, actionL = "-", noteL = paste("Lambda =", plot_output$lambda[1]), fileL = filename_dir)
    print(paste0(" - - warning: high lambda value (", plot_output$lambda[1], ")"), quote = FALSE) }
  if(plot_Manhattan & !useMan) {
    save_log(phaseL = 5L, checkL = "Manhattan plot", typeL = "insuf. data", SNPL = sum(iNA$chr | inv$chr | iNA$pos | inv$pos), allSNPs = SNPn_postQC, actionL = "-", noteL = "Insufficient known chromosome/positions to create Manhattan plot", fileL = filename_dir)
    print(" - - warning: insufficient chromosome/position values to generate Manhattan plot", quote = FALSE)
  }
  
  
  #### STEP 5: ANALYSIS complete: saving statistics data & post-QC results
  print("", quote = FALSE)
  print("Step 5: saving statistics in log file", quote = FALSE)
  flush.console()
  
  # Calculating the remaining QC statistics
  useN	<- if(SNPn_postQC > iNA_N_N + inv_N_N) TRUE else !all(iNA$N | inv$N)
  stat_N_max <- if(useN) max(dataI$N_TOTAL, na.rm = iNA_N_N+inv_N_N>0L ) else NA
  stat_N_HQ	<- if(useN & useHQ) {
    if(any(!((iNA$N | inv$N)[HQ_list]))) { max(dataI$N_TOTAL[HQ_list], na.rm = iNA_N_N+inv_N_N>0 )
    } else { NA } } else { NA }
  
  if(useFRQ & !is.na(stat_N_max)) {
    stat_Visscher <- median( 1 / (2 * dataI$EFF_ALL_FREQ * (1-dataI$EFF_ALL_FREQ) * (dataI$STDERR)^2 ), na.rm = iNA_FRQ_N+inv_FRQ_N>0L) / stat_N_HQ
    if(!is.na(stat_N_HQ) & sum(!iNA$FRQ & !inv$FRQ & HQ_list) > 9L) {
      stat_Vissc_HQ <- median( 1 / (2 * dataI$EFF_ALL_FREQ[HQ_list] * (1-dataI$EFF_ALL_FREQ[HQ_list]) * (dataI$STDERR[HQ_list])^2 ), na.rm = iNA_FRQ_N+inv_FRQ_N>0L) / stat_N_max
    } else {
      stat_Vissc_HQ <- NA
      if(is.na(stat_N_HQ))		                			save_log(4L, "Visschers stat.", "Insufficient data", SNPn_postQC, SNPn_postQC, "-", "Insufficient sample sizes to calculate Visschers statistic: too few high-quality markers", filename_dir)
      if(sum(!iNA$FRQ & !inv$FRQ & HQ_list) <= 9L)	save_log(4L, "Visschers stat.", "Insufficient data", SNPn_postQC, SNPn_postQC, "-", "Insufficient allele frequencies to calculate Visschers statistic: too few high-quality markers", filename_dir)
    }
  } else {
    if(!useFRQ) save_log(4L, "Visschers stat.", "Insufficient data", SNPn_postQC, SNPn_postQC, "-", "Insufficient allele frequencies to calculate Visschers statistic", filename_dir)
    if(is.na(stat_N_max)) save_log(4L, "Visschers stat", "Insufficient data", SNPn_postQC, SNPn_postQC, "-", "Insufficient sample sizes to calculate Visschers statistic", filename_dir)
    stat_Visscher <- NA
    stat_Vissc_HQ <- NA
  }
  
  stat_SE		<- median(dataI$STDERR)
  stat_SE_HQ		<- if(useHQ) median(dataI$STDERR[HQ_list]) else NA
  stat_skewness	<- calc_skewness(input = dataI$EFFECT)
  stat_skewness_HQ	<- if(useHQ) calc_skewness(input = dataI$EFFECT[HQ_list]) else NA
  stat_kurtosis	<- calc_kurtosis(input = dataI$EFFECT)
  stat_kurtosis_HQ	<- if(useHQ) calc_kurtosis(input = dataI$EFFECT[HQ_list]) else NA
  
  
  # Reformatting the old log entries log file
  log_old <- read.table(paste0(filename_dir, "_log.txt"), header = FALSE,
                        sep = "\t", comment.char = "", stringsAsFactors = FALSE)
  
  write.table(c(
    "*********************************",
    "",
    "\t\tQC LOG FILE",
    "",
    "*********************************",
    ""), paste0(filename_dir, "_log.txt"), append = FALSE,
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  SFL <- spreadsheet_friendly_log
  
  logCon <- file(description = paste0(filename_dir, "_log.txt"),
                 open = "a")
  
  if(SFL) {
    write.table(data.frame(
      v1 = c("Input File", "output File(s)", "QC Start Time", "QC End time", "Script version"),
      v2 = "", v3 = c(filename_input, filename, start_time, date(), "1.0-8"),
      stringsAsFactors = FALSE),
                logCon, quote = FALSE,
                sep = "\t", row.names = FALSE, col.names = FALSE)
  } else {
    write.table(format(
      data.frame(
        v1 = c("Input File", "output File(s)", "QC Start Time", "QC End time", "Script version"),
        v2 = "\t: ", v3 = c(filename_input, filename, start_time, date(), "1.0-8"),
        stringsAsFactors = FALSE), justify = "left"),
                logCon, quote = FALSE,
                sep = "", row.names = FALSE, col.names = FALSE)
  }
  
  write.table(c("", "",
                "****************************************",
                "\t1. Log entries generated during QC",
                "****************************************",
                ""), logCon, quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  write.table(if(SFL) log_old else format(log_old, justify = "left"),
              logCon, quote = FALSE,
              sep = if(SFL) "\t" else "   ", row.names = FALSE, col.names = FALSE)
  rm(log_old)
  
  write.table(c("", "",
                "*******************************************",
                "\t2. Summary statistics of removed SNPs",
                "*******************************************",
                ""), logCon, quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  
  
  del_table <- data.frame(
    # Statistics of remove 1: monomorphic, Y- & M-chr SNPs
    name1 = c("SNPs in uploaded data", "Removed", "> Monomorphic", ">> missing non-effect allele", ">> invalid non-effect allele", ">> allele frequency = 0", ">> identical alleles*", if(remove_X) "> X chromosome" else ">", if(remove_Y) "> Y chromosome" else ">", if(remove_XY) "> XY chromosome" else ">", if(remove_M) "> M chromosome" else ">", "", "", ""),
    N1 = c(SNPn_input, remove_L0_N, monomorp_N, iNA_al2_N, inv_al2_N, low_FRQ_N, same_al_N, chr_X_N, chr_Y_N, chr_XY_N, chr_M_N, NA, NA, NA),
    perc1 = numeric(length = 14),
    em4 = character(length = 14),
    # Statistics of remove 2: missing & invalid crucial variables
    name2 = c("SNPs in initial QC", "Removed", "> Marker name", ">> missing", ">> duplicate", "> Effect allele", ">> missing", ">> invalid", "> Effect size (missing)", ">> low imputation quality", "> Standard error", ">> missing", ">>> low imputation quality", ">> invalid"),
    N2 = c(SNPn_preQC, remove_L1_N, iNA_mar_N + SNPn_preQC_dupli, iNA_mar_N, SNPn_preQC_dupli, iNA_al1_N + inv_al1_N, iNA_al1_N, inv_al1_N, iNA_eff_N, iNA_eff_poor, iNA_se_N + inv_se_N, iNA_se_N, iNA_se_poor,  inv_se_N),
    perc2 = numeric(length = 14),
    em8 = character(length = 14),
    # Statistics of remove 3: allele mismatches & FINAL statistics
    name3 = c(character(length = 10), "Final dataset", "SNPs", "> genotyped", "> imputed"),
    N3 = c(rep(NA, 11), SNPn_postQC, SNPn_postQC_geno, SNPn_postQC_imp),
    perc3 = numeric(length = 14),
    stringsAsFactors = FALSE)
  del_table$perc1 <- round(100 * del_table$N1 / SNPn_input, digits = 2)
  del_table$perc2 <- round(100 * del_table$N2 / SNPn_preQC, digits = 2)
  del_table$perc3 <- round(100 * del_table$N3 / SNPn_postQC, digits = 2)
  
  if(remove_mismatches_std | remove_mismatches_alt | remove_diffEAF_std | remove_diffEAF_alt) {
    del_table$name3[1:2] <- c("SNPs in allele check", "Removed")
    del_table$N3[1] <- SNPn_ref
    del_table$N3[2] <- remove_L2_N
    
    if(remove_mismatches_std & use_allele_std){
      del_table$name3[3] <- paste("> mismatch to", allele_name_std)
      del_table$N3[3] <- allele_out_std$n_mismatch    
    } else { del_table$name3[3] <- ">" }
    if(remove_mismatches_alt & use_allele_alt){
      del_table$name3[4] <- paste("> mismatch to", allele_name_alt)
      del_table$N3[4] <- allele_out_alt$n_mismatch    
    } else { del_table$name3[4] <- ">" }
    
    if(remove_diffEAF_std & use_allele_std){
      del_table$name3[5] <- paste("> deviating FRQ from", allele_name_std)
      del_table$N3[5] <- allele_out_std$n_diffEAF    
    } else { del_table$name3[5] <- ">" }
    if(remove_diffEAF_alt & use_allele_alt){
      del_table$name3[6] <- paste("> deviating FRQ from", allele_name_alt)
      del_table$N3[6] <- allele_out_alt$n_diffEAF    
    } else { del_table$name3[6] <- ">" }
    
    del_table$perc3[1:6] <- round(100 * del_table$N3[1:6] / SNPn_ref, digits = 2)
  }
  
  # adding the header line
  del_table <- rbind(character(length = 11), del_table)
  del_table[1, ] <- if(remove_mismatches_std | remove_mismatches_alt | remove_diffEAF_std | remove_diffEAF_alt) {
    c("", "N", "%", "", "", "N", "%", "", "", "N", "%")
  } else { c("", "N", "%", "", "", "N", "%", "", "", "", "") }
  del_table[12, c("N3", "perc3")] <- c("N", "%")
  
  # Changing NA's to empty cells
  del_table$N1 <- ifelse(is.na(del_table$N1), "", del_table$N1)
  del_table$perc1 <- ifelse(is.na(del_table$perc1), "", del_table$perc1)
  del_table$N3 <- ifelse(is.na(del_table$N3), "", del_table$N3)
  del_table$perc3 <- ifelse(is.na(del_table$perc3), "", del_table$perc3)
  
  write.table(if(SFL) del_table else format(del_table, justify = "left"),
              logCon, quote = FALSE,
              sep = if(SFL) "\t" else "   ", row.names = FALSE, col.names = FALSE)
  rm(del_table)
  
  
  write.table(c("* Not including SNPs with allele frequency = 0",
                "     NB: monomorphic & Y- & M-chromosome SNPs are removed before the analysis starts.",
                "     The pre-QC values in the tables below refer to the dataset after these SNPs have been",
                "     removed, but before further exclusions have taken place.",
                "", "",
                "********************************************************************",
                "\t3. Summary statistics before and after quality check procedure",
                "********************************************************************",
                ""),
              logCon, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  stat_table <- data.frame(variable = character(18),
                           NA_pre = character(18), inv_pre = character(18), unus_pre = character(18),
                           NA_post = character(18), inv_post = character(18), unus_post = character(18),
                           Qmin = character(18), Q25 = character(18), Qmean = character(18), Qmedian = character(18), Q75 = character(18), Qmax = character(18),
                           stringsAsFactors = FALSE)
  
  stat_table[1, 2] <- "Pre-QC"
  stat_table[1, 5] <- "Post-QC"
  stat_table[2,  ] <- c("Statistics", "missing", "invalid", "%unusable", "missing", "invalid", "%unusable", "min", "25%", "mean", "median", "75%", "max")
  
  stat_save <- function(stat_name, stat_col, old_N, old_NA, old_inv,
                        new_N = length(stat_col), new_NA = sum(is.na(stat_col)) - new_inv, new_inv,
                        no_quantiles = FALSE, round_min = FALSE) {
    if(no_quantiles | new_NA + new_inv == new_N) {
      return(c(stat_name, old_NA, old_inv, round(100*(old_NA + old_inv)/old_N, digits = 2), new_NA, new_inv, round(100*(new_NA + new_inv)/new_N, digits = 2), character(6) ))
    } else {
      quant <- quantile(stat_col, na.rm = new_NA + new_inv > 0L, names = FALSE)
      return(c(stat_name, old_NA, old_inv, round(100*(old_NA + old_inv)/old_N, digits = 2), new_NA, new_inv, round(100*(new_NA + new_inv)/new_N, digits = 2), if(round_min) round(quant[1], digits = 2) else signif(quant[1], digits = 2), round(c(mean(stat_col, na.rm = new_NA + new_inv > 0L), quant)[c(3,1,4,5,6)], digits = 4)))
    }
  }
  stat_select_save <- function(stat_name, select_name,
                               stat_col, NA_col = is.na(stat_col) & !inv_col, inv_col, select_col,
                               old_N, old_NA, old_inv,
                               new_N = sum(select_col), new_NA = sum(select_col & NA_col), new_inv = sum(select_col & inv_col),
                               round_min = FALSE) {
    if(old_N == 0L) {
      return(c(paste(stat_name, "-", select_name), "-", "-", "-", "-", "-", "-", character(6)))
    } else {
      if(new_N == 0L) {
        return(c(paste(stat_name, "-", select_name), old_NA, old_inv, round(100*(old_NA + old_inv)/old_N, digits = 2), "-", "-", "-", character(6)))
      } else {
        new_inv <- sum(select_col & inv_col)
        new_NA  <- sum(select_col & NA_col)
        if(new_NA + new_inv == new_N) {
          return(c(paste(stat_name, "-", select_name), old_NA, old_inv, round(100*(old_NA + old_inv)/old_N, digits = 2), new_NA, new_inv, round(100*(new_NA + new_inv)/new_N, digits = 2), character(6)))
        } else {
          quant <- quantile(stat_col[select_col], na.rm = new_NA + new_inv > 0L, names = FALSE)
          return(c(paste(stat_name, "-", select_name), old_NA, old_inv, round(100*(old_NA + old_inv)/old_N, digits = 2), new_NA, new_inv, round(100*(new_NA + new_inv)/new_N, digits = 2), if(round_min) round(quant[1], digits = 2) else signif(quant[1], digits = 2), round(c(mean(stat_col[select_col], na.rm = new_NA + new_inv > 0L), quant)[c(3,1,4,5,6)], digits = 4)))
        } } }
  }
  
  stat_table[ 3, ] <- stat_save(stat_name = "effect size", stat_col = dataI$EFFECT, old_N = SNPn_preQC, old_NA = iNA_eff_N, old_inv = 0, new_N = SNPn_postQC, new_NA = 0, new_inv = 0)
  stat_table[ 4, ] <- stat_save(stat_name = "SE", stat_col = dataI$STDERR, old_N = SNPn_preQC, old_NA = iNA_se_N, old_inv = inv_se_N, new_N = SNPn_postQC, new_NA = 0, new_inv = 0)
  stat_table[ 5, ] <- stat_save(stat_name = "p-value", stat_col = dataI$PVALUE, old_N = SNPn_preQC, old_NA = iNA_p_N, old_inv = inv_p_N, new_N = SNPn_postQC, new_NA = sum(iNA$p), new_inv = sum(inv$p))
  stat_table[ 6, ] <- stat_save(stat_name = "allele frequency", stat_col = dataI$EFF_ALL_FREQ, old_N = SNPn_preQC, old_NA = iNA_FRQ_N, old_inv = inv_FRQ_N, new_N = SNPn_postQC, new_NA = sum(iNA$FRQ), new_inv = sum(inv$FRQ))
  stat_table[ 7, ] <- stat_save(stat_name = "HWE p-value", stat_col = dataI$HWE_PVAL, old_N = SNPn_preQC, old_NA = iNA_HWE_N, old_inv = inv_HWE_N, new_N = SNPn_postQC, new_NA = sum(iNA$HWE), new_inv = sum(inv$HWE))
  stat_table[ 8, ] <- stat_select_save(stat_name = "HWE p-value", select_name = "genotyped", stat_col = dataI$HWE_PVAL, NA_col = iNA$HWE, inv_col = inv$HWE,
                                       select_col = geno_list,  old_N = SNPn_preQC_geno, old_NA = iNA_HWE_g, old_inv = inv_HWE_g, new_N = SNPn_postQC_geno)
  stat_table[ 9, ] <- stat_select_save(stat_name = "HWE p-value", select_name = "imputed", stat_col = dataI$HWE_PVAL, NA_col = iNA$HWE, inv_col = inv$HWE,
                                       select_col =  imp_list,  old_N = SNPn_preQC_imp, old_NA = iNA_HWE_i, old_inv = inv_HWE_i, new_N = SNPn_postQC_imp)
  stat_table[10, ] <- stat_save(stat_name = "Call rate", stat_col = dataI$CALLRATE, old_N = SNPn_preQC, old_NA = iNA_cal_N, old_inv = inv_cal_N, new_N = SNPn_postQC, new_NA = sum(iNA$cal), new_inv = sum(inv$cal))
  stat_table[11, ] <- stat_select_save(stat_name = "Call rate", select_name = "genotyped", stat_col = dataI$CALLRATE, NA_col = iNA$cal, inv_col = inv$cal,
                                       select_col = geno_list,	old_N = SNPn_preQC_geno, old_NA = iNA_cal_g, old_inv = inv_cal_g, new_N = SNPn_postQC_geno)
  stat_table[12, ] <- stat_select_save(stat_name = "Call rate", select_name = "imputed", stat_col = dataI$CALLRATE, NA_col = iNA$cal, inv_col = inv$cal,
                                       select_col =	imp_list,	old_N = SNPn_preQC_imp, old_NA = iNA_cal_i, old_inv = inv_cal_i, new_N = SNPn_postQC_imp)
  stat_table[13, ] <- stat_save(stat_name = "Sample size", stat_col = dataI$N_TOTAL, old_N = SNPn_preQC, old_NA = iNA_N_N, old_inv = inv_N_N, new_N = SNPn_postQC, new_NA = sum(iNA$N), new_inv = sum(inv$N), round_min = TRUE)
  stat_table[14, ] <- stat_select_save(stat_name = "Sample size", select_name = "genotyped", stat_col = dataI$N_TOTAL, NA_col = iNA$N, inv_col = inv$N,
                                       select_col = geno_list,	old_N = SNPn_preQC_geno, old_NA = iNA_N_g, old_inv = inv_N_g, new_N = SNPn_postQC_geno, round_min = TRUE)
  stat_table[15, ] <- stat_select_save(stat_name = "Sample size", select_name = "imputed", stat_col = dataI$N_TOTAL, NA_col = iNA$N, inv_col = inv$N,
                                       select_col =	imp_list,	old_N = SNPn_preQC_imp, old_NA = iNA_N_i, old_inv = inv_N_i, new_N = SNPn_postQC_imp, round_min = TRUE)
  stat_table[16, ] <- stat_save(stat_name = "Imputation quality", stat_col = dataI$IMP_QUALITY, old_N = SNPn_preQC, old_NA = iNA_impQ_N, old_inv = inv_impQ_N, new_N = SNPn_postQC, new_NA = sum(iNA$impQ), new_inv = sum(inv$impQ))
  stat_table[17, ] <- stat_select_save(stat_name = "Imp. quality", select_name = "genotyped", stat_col = dataI$IMP_QUALITY, NA_col = iNA$impQ, inv_col = inv$impQ,
                                       select_col = geno_list,	old_N = SNPn_preQC_geno, old_NA = iNA_impQ_g, old_inv = inv_impQ_g, new_N = SNPn_postQC_geno)
  stat_table[18, ] <- stat_select_save(stat_name = "Imp. quality", select_name = "imputed", stat_col = dataI$IMP_QUALITY, NA_col = iNA$impQ, inv_col = inv$impQ,
                                       select_col =	imp_list,	old_N = SNPn_preQC_imp, old_NA = iNA_impQ_i, old_inv = inv_impQ_i, new_N = SNPn_postQC_imp)
  write.table(if(SFL) stat_table else format(stat_table, justify = "left"),
              logCon, quote = FALSE,
              sep = if(SFL) "\t" else "   ", row.names = FALSE, col.names = FALSE)
  rm(stat_table)
  
  write.table(c("", "",
                "************************",
                "\t4. Allele matching",
                "************************",
                ""),
              logCon, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  all_table0 <- data.frame(
    names = c("SNPs", "> negative strand", "> strand switch (SS)", ">> double SS",
              ">> MISMATCH", ">>> mismatch & double SS",
              "> flipped", "> ambiguous", ">> MAF between 35 & 65%", "> deviating allele frequency"),
    N = integer(length = 10), perc = numeric(length = 10), stringsAsFactors = FALSE)
  
  if(use_allele_std & SNPn_ref_std > 0L) {
    all_table1 <- all_table0
    all_table1$N <- c(allele_out_std$n_SNPs, if(strand_ref$minus == 0L) 0L else allele_out_std$n_negative_strand, allele_out_std$n_strandswitch, if(strand_ref$minus == 0L) 0L else allele_out_std$n_negative_switch,
                      allele_out_std$n_mismatch, if(strand_ref$minus == 0L) 0L else allele_out_std$n_negative_mismatch,
                      allele_out_std$n_flipped, allele_out_std$n_ambiguous, allele_out_std$n_suspect, allele_out_std$n_diffEAF)
    all_table1$perc <- round(100 * all_table1$N/SNPn_ref, digits = 2)
    all_table0$N <- all_table1$N
    
    all_table1$N    <- ifelse(is.na(all_table1$N   ), "-", as.character(all_table1$N   ))
    all_table1$perc <- ifelse(is.na(all_table1$perc), "-", as.character(all_table1$perc))
  } else { all_table1 <- cbind(all_table0$names, matrix(data = "-", nrow = 10, ncol = 2)) }
  
  if(use_allele_alt & SNPn_ref_alt > 0L) {
    all_table2 <- all_table0
    all_table2$N <- c(allele_out_alt$n_SNPs, if(strand_ref$minus == 0L) 0L else allele_out_alt$n_negative_strand, allele_out_alt$n_strandswitch, if(strand_ref$minus == 0L) 0L else allele_out_alt$n_negative_switch,
                      allele_out_alt$n_mismatch, if(strand_ref$minus == 0L) 0L else allele_out_alt$n_negative_mismatch,
                      allele_out_alt$n_flipped, allele_out_alt$n_ambiguous, allele_out_alt$n_suspect, allele_out_alt$n_diffEAF)
    all_table2$perc <- round(100 * all_table2$N/SNPn_ref, digits = 2)
    all_table0$N <- all_table0$N + all_table2$N
    
    all_table2$N    <- ifelse(is.na(all_table2$N   ), "-", as.character(all_table2$N   ))
    all_table2$perc <- ifelse(is.na(all_table2$perc), "-", as.character(all_table2$perc))
  } else { all_table2 <- cbind(all_table0$names, matrix(data = "-", nrow = 10, ncol = 2)) }
  
  if(SNPn_ref_new > 0L) {
    all_table3 <- all_table0
    all_table3$N <- c(SNPn_ref_new, if(strand_ref$minus == 0L) 0L else allele_out_new$n_negative_strand, NA, NA,
                      NA, NA, allele_out_new$n_flipped, allele_out_new$n_ambiguous, NA, NA)
    all_table3$perc <- round(100 * all_table3$N/SNPn_ref, digits = 2)
    all_table0$N[c(1,2,7,8)] <- all_table0$N[c(1,2,7,8)] + all_table3$N[c(1,2,7,8)]
    
    all_table3$N    <- ifelse(is.na(all_table3$N   ), "-", as.character(all_table3$N   ))
    all_table3$perc <- ifelse(is.na(all_table3$perc), "-", as.character(all_table3$perc))
  } else { all_table3 <- cbind(all_table0$names, c("0", rep("-", 9)), c("0", rep("-", 9))) }
  
  all_table0$perc <- round(100 * all_table0$N/SNPn_ref, digits = 2)
  all_table0$N    <- ifelse(is.na(all_table0$N   ), "-", as.character(all_table0$N   ))
  all_table0$perc <- ifelse(is.na(all_table0$perc), "-", as.character(all_table0$perc))
  
  
  # Combining tables & adding header line
  write.table(if(SFL) {
    rbind(c("Combined", "N", "%", "", allele_name_std, "N", "%", "", allele_name_alt, "N", "%", "", if(update_alt) paste("Updated", allele_name_alt) else "Other SNPs", "N", "%"),
          cbind(all_table0, character(length = 10), all_table1, character(length = 10), all_table2, character(length = 10), all_table3, stringsAsFactors = F))
  } else {
    format(rbind(c("Combined", "N", "%", "", allele_name_std, "N", "%", "", allele_name_alt, "N", "%", "", if(update_alt) paste("Updated", allele_name_alt) else "Other SNPs", "N", "%"),
                 cbind(all_table0, character(length = 10), all_table1, character(length = 10), all_table2, character(length = 10), all_table3, stringsAsFactors = F)), justify = "left") },
              logCon, quote = FALSE, sep = if(SFL) "\t" else "   ", row.names = FALSE, col.names = FALSE)
  rm(all_table0, all_table1, all_table2, all_table3)
  
  
  write.table(c("", "",
                "**********************",
                "\t5. QC statistics",
                "**********************",
                ""),
              logCon, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  
  stat_table <- data.frame(allname = c(allele_name_std, "> ambiguous SNPs", "> non-ambiguous SNPs", allele_name_alt, "> ambiguous SNPs", "> non-ambiguous SNPs"),
                           allcor = "-",
                           em1 = character(length = 6),
                           name1 = c("r", "Lambda", "> genotyped", "> imputed", "", ""),
                           stat1 = c(round(outcome_P, digits = 3), round(plot_output$lambda, digits = 3), "", ""),
                           em2 = character(length = 6),
                           name2 = c("SE median", "Skewness", "Kurtosis", "Visscher's stat.", "* high-quality SNPs", "   only"),
                           stat2all = c(signif(stat_SE, digits = 3), round(c(stat_skewness, stat_kurtosis, stat_Visscher), digits = 3), "", ""),
                           stat2_HQ = c(signif(stat_SE_HQ, digits = 3), round(c(stat_skewness_HQ, stat_kurtosis_HQ, stat_Vissc_HQ), digits = 3), "", ""),
                           em3 = character(length = 6),
                           Nname = c("Negative strand SNPs", "High-quality SNPs", "Corrected p-values", "Extreme p-values", "", ""),
                           NN = c(strand_post$minus, SNPn_postQC_HQ, calc_p_newN, low_p_newN, NA, NA),
                           Nperc = numeric(length = 6), stringsAsFactors = FALSE)
  stat_table$Nperc <- round(100*stat_table$NN/SNPn_postQC, digits = 2)		
  stat_table[5:6, c("NN", "Nperc")] <- "" 
  if(!calculate_missing_p) { stat_table[3:4, c("NN", "Nperc")] <- "-" }
  if(use_allele_std | use_allele_alt) {
    stat_table$allcor <- as.character(round(c(allele_out_std$FRQ_cor, allele_out_std$FRQ_cor_ambiguous, allele_out_std$FRQ_cor_nonambi, allele_out_alt$FRQ_cor, allele_out_alt$FRQ_cor_ambiguous, allele_out_alt$FRQ_cor_nonambi), digits = 3))
    if(!check_ambiguous_alleles) { stat_table$allcor[c(2,3,5,6)] <- "-" }
    if(!use_allele_alt) { stat_table$allcor[c(4,5,6)] <- "-" }
  }
  
  write.table(if(SFL) {
    rbind(c("Allele-frequency correlation", "", "", "P-value corr.", "", "", "QC stats", "all SNPs", "*HQ SNPs", "", "Other", "N", "%"), stat_table)
  } else {
    format(rbind(c("Allele-frequency correlation", "", "", "P-value corr.", "", "", "QC stats", "all SNPs", "*HQ SNPs", "", "Other", "N", "%"),
                 stat_table), justify = "left") },
              logCon, quote = FALSE, sep = if(SFL) "\t" else "   ", row.names = FALSE, col.names = FALSE)
  rm(stat_table)
  
  if(ignore_impstatus) {
    write.table(c("", "",
                  "************************",
                  "\t6. QQ plot filters",
                  "************************",
                  "",
                  "NOTE - ignore_impstatus was TRUE: HWE-p, callrate and imputation-Q filters (HQ & QQ) were applied to all SNPs",
                  "",
                  ""),
                logCon, quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(c("", "",
                  "************************",
                  "\t6. QQ plot filters",
                  "************************",
                  "",
                  "NOTE - Ignore_impstatus was FALSE: HWE-p and callrate filters (HQ & QQ) were applied to genotyped SNPs only;",
                  "       imputation-Q filter to imputed SNPs only.",
                  ""),
                logCon, quote = FALSE, row.names = FALSE, col.names = FALSE) }
  
  
  FRQ_table <- data.frame(names = character(length = 5), N = character(5), perc = character(5), em = character(length = 5), stringsAsFactors = FALSE)
  if(is.null(plot_output$FRQfilter_N)) { FRQ_table$names[1] <- "no filter applied"
  } else {
    FRQ_table$names[1:length(plot_output$FRQfilter_N)] <- plot_output$FRQfilter_names
    FRQ_table$N[1:length(plot_output$FRQfilter_N)] <- as.character(plot_output$FRQfilter_N)
    FRQ_table$perc[1:length(plot_output$FRQfilter_N)] <- as.character(round(100*plot_output$FRQfilter_N/SNPn_postQC, digits = 2))
  }
  
  HWE_table <- data.frame(names = character(length = 5), N = character(5), perc = character(5), em = character(length = 5), stringsAsFactors = FALSE)
  if(is.null(plot_output$HWEfilter_N)) { HWE_table$names[1] <- "no filter applied"
  } else {
    HWE_table$names[1:length(plot_output$HWEfilter_N)] <- plot_output$HWEfilter_names
    HWE_table$N[1:length(plot_output$HWEfilter_N)] <- as.character(plot_output$HWEfilter_N)
    HWE_table$perc[1:length(plot_output$HWEfilter_N)] <- as.character(round(100*plot_output$HWEfilter_N/SNPn_postQC, digits = 2))
  }
  
  cal_table <- data.frame(names = character(length = 5), N = character(5), perc = character(5), em = character(length = 5), stringsAsFactors = FALSE)
  if(is.null(plot_output$calfilter_N)) { cal_table$names[1] <- "no filter applied"
  } else {
    cal_table$names[1:length(plot_output$calfilter_N)] <- plot_output$calfilter_names
    cal_table$N[1:length(plot_output$calfilter_N)] <- as.character(plot_output$calfilter_N)
    cal_table$perc[1:length(plot_output$calfilter_N)] <- as.character(round(100*plot_output$calfilter_N/SNPn_postQC, digits = 2))
  }
  
  imp_table <- data.frame(names = character(length = 5), N = character(5), perc = character(5), stringsAsFactors = FALSE)
  if(is.null(plot_output$impfilter_N)) { imp_table$names[1] <- "no filter applied"
  } else {
    imp_table$names[1:length(plot_output$impfilter_N)] <- plot_output$impfilter_names
    imp_table$N[1:length(plot_output$impfilter_N)] <- as.character(plot_output$impfilter_N)
    imp_table$perc[1:length(plot_output$impfilter_N)] <- as.character(round(100*plot_output$impfilter_N/SNPn_postQC, digits = 2))
  }
  
  write.table(if(SFL) {
    rbind(c("Allele frequency", "N", "%", "", "HWE p-value", "N", "%", "", "Call rate", "N", "%", "", "Imp. quality", "N", "%"),
          cbind(FRQ_table, HWE_table, cal_table, imp_table))  
  } else {
    format(rbind(c("Allele frequency", "N", "%", "", "HWE p-value", "N", "%", "", "Call rate", "N", "%", "", "Imp. quality", "N", "%"),
                 cbind(FRQ_table, HWE_table, cal_table, imp_table)), justify = "left") },
              logCon, quote = FALSE, sep = if(SFL) "\t" else "   ", row.names = FALSE, col.names = FALSE)
  rm(FRQ_table, HWE_table, cal_table, imp_table)
  
  
  write.table(c("", "",
                "*********************************",
                "\t7. Chromosome & Allele data",
                "*********************************",
                ""),
              logCon, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  ### Tables counting chromosomes & alleles
  chr_table <- data.frame(Chr = c("NA", "Invalid", 1:22, "23 ( X)",
                                  "24 ( Y)", "25 (XY)", "26 ( M)"),
                          N1 = integer(length = 28),
                          perc1= numeric(length = 28),
                          empty= character(length = 28),
                          all = character(length = 28),
                          N2 = integer(length = 28),
                          perc2= numeric(length = 28),
                          stringsAsFactors = FALSE)
  
  chr_table$N1[1] <- sum(iNA$chr)
  chr_table$N1[2] <- sum(inv$chr)
  na_rm_chr <- chr_table$N1[1] + chr_table$N1[2] > 0L
  for(ci in 1:22) { chr_table$N1[ci + 2L] <- sum(dataI$CHR == ci, na.rm = na_rm_ch) }
  chr_table$N1[25] <- if(remove_X ) NA else sum(dataI$CHR == 23L, na.rm = na_rm_ch)
  chr_table$N1[26] <- if(remove_Y ) NA else sum(dataI$CHR == 24L, na.rm = na_rm_ch)
  chr_table$N1[27] <- if(remove_XY) NA else sum(dataI$CHR == 25L, na.rm = na_rm_ch)
  chr_table$N1[28] <- if(remove_M ) NA else sum(dataI$CHR == 26L, na.rm = na_rm_ch)
  chr_table$perc1 <- round(100*chr_table$N1/SNPn_postQC, digits = 2)
  
  if(remove_X) {
    chr_table$N1[25] <- "removed"
    chr_table$perc1[25] <- paste0("(", chr_X_N, " SNPs)") }
  if(remove_Y) {
    chr_table$N1[26] <- "removed"
    chr_table$perc1[26] <- paste0("(", chr_Y_N, " SNPs)") }
  if(remove_XY) {
    chr_table$N1[27] <- "removed"
    chr_table$perc1[27] <- paste0("(",chr_XY_N, " SNPs)") }
  if(remove_M) {
    chr_table$N1[28] <- "removed"
    chr_table$perc1[28] <- paste0("(", chr_M_N, " SNPs)") }
  
  chr_table$all[1:10] <- c("A", "T", "C", "G", "", "Other allele", "A", "T", "C", "G")
  chr_table$N2[1] <- sum(dataI$EFFECT_ALL == "A")
  chr_table$N2[2] <- sum(dataI$EFFECT_ALL == "T")
  chr_table$N2[3] <- sum(dataI$EFFECT_ALL == "C")
  chr_table$N2[4] <- sum(dataI$EFFECT_ALL == "G")
  chr_table$N2[7] <- sum(dataI$OTHER_ALL == "A")
  chr_table$N2[8] <- sum(dataI$OTHER_ALL == "T")
  chr_table$N2[9] <- sum(dataI$OTHER_ALL == "C")
  chr_table$N2[10]<- sum(dataI$OTHER_ALL == "G")
  chr_table$perc2[1:10] <- round(100*chr_table$N2[1:10]/SNPn_postQC, digits = 2)
  chr_table[c(5, 6, 11:28), c(6,7)] <- ""
  
  write.table(if(SFL) {
    rbind(c("Chromosome", "N", "%", "", "Effect allele", "N", "%"), chr_table)
  } else {
    format(rbind(c("Chromosome", "N", "%", "", "Effect allele", "N", "%"), chr_table), justify = "left") },
              logCon, quote = FALSE, sep = if(SFL) "\t" else "   ", row.names = FALSE, col.names = FALSE)
  rm(chr_table)
  
  
  write.table(c("", "",
                "*********************************",
                "\t8. Missing & invalid values",
                "*********************************",
                ""),
              logCon, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  error_table1 <- data.frame(names = c("Duplicate SNPs", "Imputation status", "Strand", "> missing", "> invalid", "Chromosome", "> missing", "> invalid", "Position", "> missing", "> invalid", "Invalid SNPs"),
                             N = c(SNPn_preQC_dupli, SNPn_preQC	- (SNPn_preQC_geno+ SNPn_preQC_imp), strand_pre$missing + strand_pre$invalid, strand_pre$missing, strand_pre$invalid, iNA_chr_N + inv_chr_N, iNA_chr_N, inv_chr_N, iNA_pos_N + inv_pos_N, iNA_pos_N, inv_pos_N, SNPn_preQC_inv),
                             perc = numeric(length = 12), stringsAsFactors = FALSE)
  error_table1$perc <- round(100*error_table1$N / SNPn_preQC, digits = 2)
  error_table2 <- error_table1
  error_table2[1,1] <- "Extreme p-values"
  error_table2$N <- c(low_p_newN, SNPn_postQC - (SNPn_postQC_geno+ SNPn_postQC_imp), strand_post$missing + strand_post$invalid, strand_post$missing, strand_post$invalid, sum(iNA$chr | inv$chr), sum(iNA$chr), sum(inv$chr), sum(iNA$pos | inv$pos), sum(iNA$pos), sum(inv$pos), SNPn_postQC_inv)
  error_table2$perc <- round(100*error_table2$N / SNPn_postQC, digits = 2)
  
  write.table(if(SFL) {
    rbind(c("Pre-QC", "N", "%", "", "Post-QC", "N", "%"),
          cbind(error_table1, character(length = 12), error_table2))
  } else {
    format(rbind(c("Pre-QC", "N", "%", "", "Post-QC", "N", "%"),
                 cbind(error_table1, character(length = 12), error_table2)),
           justify = "left") },
              logCon, quote = FALSE, sep = if(SFL) "\t" else "   ", row.names = FALSE, col.names = FALSE)
  write.table("   * Note: SNPs whose imputation status is missing (i.e. not invalid) are also counted as invalid.",
              logCon, col.names=FALSE, row.names=FALSE, quote=FALSE)
  rm(error_table1, error_table2)
  
  
  write.table(c("", "",
                "***************************",
                "\t9. Settings of the QC",
                "***************************",
                ""),
              logCon, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  settings_table <- data.frame(
    col1 = c("input file", "output file", "dir data", "dir output", "dir refs", ""),
    col2 = c(filename_input, filename, dir_data, dir_output, dir_references, ""),
    stringsAsFactors = FALSE)
  
  write.table(if(SFL) settings_table else format(settings_table, justify = "left"),
              logCon, col.names=FALSE, row.names=FALSE, quote=FALSE, sep=if(SFL) "\t" else "   ")
  rm(settings_table)
  # Table already includes an additional row to create whiteline
  
  filter_table <- data.frame(
    chr_names = c(" X chrom.", " Y chrom.", "XY chrom.", " M chrom."),
    chr_val = c(remove_X, remove_Y, remove_XY, remove_M),
    em1 = "   ",
    fil_names = c("Allele FRQ", "HWE p-value", "Callrate", "Imp. Quality"),
    threshold = c(useFRQ_threshold, useHWE_threshold, useCal_threshold, useImp_threshold),
    fNA = logical(length = 4),
    HQ = c(NA, NA, NA, NA),
    QQ = c("-", "-", "-", "-"),
    stringsAsFactors = FALSE)
  if(useFRQ) {
    if(!is.null(c(HQfilter_FRQ, QQfilter_FRQ)))
      filter_table[1,6] <- NAfilter_FRQ
    if(!is.null(HQfilter_FRQ)) filter_table[1,7] <- HQfilter_FRQ
    if(!is.null(QQfilter_FRQ)) filter_table[1,8] <- paste(QQfilter_FRQ, collapse = ", ")
  } else { filter_table[1,8] <- "Insufficient data to meet threshold" }
  if(useHWE) {
    if(!is.null(c(HQfilter_HWE, QQfilter_HWE)))
      filter_table[2,6] <- NAfilter_HWE
    if(!is.null(HQfilter_HWE)) filter_table[2,7] <- HQfilter_HWE
    if(!is.null(QQfilter_HWE)) filter_table[2,8] <- paste(QQfilter_HWE, collapse = ", ")
  } else { filter_table[2,8] <- "Insufficient data to meet threshold" }
  if(useCal) {
    if(!is.null(c(HQfilter_cal, QQfilter_cal)))
      filter_table[3,6] <- NAfilter_cal
    if(!is.null(HQfilter_cal)) filter_table[3,7] <- HQfilter_cal
    if(!is.null(QQfilter_cal)) filter_table[3,8] <- paste(QQfilter_cal, collapse = ", ")
  } else { filter_table[3,8] <- "Insufficient data to meet threshold" }
  if(useImp) {
    if(!is.null(c(HQfilter_imp, QQfilter_imp)))
      filter_table[4,6] <- NAfilter_imp
    if(!is.null(HQfilter_imp)) filter_table[4,7] <- HQfilter_imp
    if(!is.null(QQfilter_imp)) filter_table[4,8] <- paste(QQfilter_imp, collapse = ", ")
  } else { filter_table[4,8] <- "Insufficient data to meet threshold" }
  
  write.table(if(SFL) {
    rbind(c("> filters", "", "", "", "threshold*", "NA", "HQ", "QQ"), filter_table)
  } else {
    format(rbind(c("> filters", "", "", "", "threshold*", "NA", "HQ", "QQ"),
                 filter_table), justify = "left") },
              logCon, quote = FALSE, sep = if(SFL) "\t" else "   ", row.names = FALSE, col.names = FALSE)
  write.table(c("* the value used in phase 4 (opt. after multiplication with the remaining number of SNPs)", ""),
              logCon, quote = FALSE, row.names = FALSE, col.names = FALSE)
  rm(filter_table)
  
  settings_table <- rbind(c("> other settings", "", "", "", "", "", "> plots", "", "", "", "", "", "> allele references", "", "", "", ""),
                          data.frame(
                            c1 = c("remove_mismatches_std", "remove_mismatches_alt", "remove_diffEAF_std", "remove_diffEAF_alt", "threshold_diffEAF", "check_ambiguous_alleles", "return_HQ_effectsizes", "calculate_missing_p"),
                            c2 = c(remove_mismatches_std, remove_mismatches_alt, remove_diffEAF_std, remove_diffEAF_alt, as.character(threshold_diffEAF), check_ambiguous_alleles, return_HQ_effectsizes, calculate_missing_p),
                            c3 = "   ",
                            c4 = c("ignore impS", "min imp", "max imp", "imputed T", "imputed F", "imputed NA","",""),
                            c5 = c(ignore_impstatus,  as.character(minimal_impQ_value), as.character(maximal_impQ_value), paste(imputed_T, collapse = ", "), paste(imputed_F, collapse = ", "), paste(imputed_NA, collapse = ", "),"",""),
                            c6 = "   ",
                            c7 = c("make plots", "only plot if threshold", "histograms", "QQ plot", "Manhattan", "intensity","",""),
                            c8 = c(make_plots, only_plot_if_threshold, plot_histograms, plot_QQ, plot_Manhattan, plot_intensity,"",""),
                            c9 = "   ",
                            c10= c("threshold FRQ test" , "threshold p-test", "QQ bands", "plot p cutoff", "*Manhattan threshold", "","",""),
                            c11= c(threshold_allele_freq_correlation, threshold_p_correlation, as.character(plot_QQ_bands), plot_cutoff_p, useMan_threshold, "","",""),
                            c12= "   ",
                            c13= c("update alt", "update savename", "save as rdata", "backup alt", "", "","",""),
                            c14= if(update_alt) { c(TRUE, update_savename, update_as_rdata, backup_alt, "", "","","")
                            } else { c(FALSE, "-", "-", "-", "", "", "", "" ) },
                            c15= "   ",
                            c16= c("allele ref std input", "allele ref std name", "allele ref alt input", "allele ref alt name", "", "","",""),
                            c17= c(settings_allele_std, allele_name_std, settings_allele_alt, allele_name_alt, "", "","",""),
                            stringsAsFactors = FALSE))
  write.table(if(SFL) settings_table else format(settings_table, justify = "left"),
              logCon, quote = FALSE, sep = if(SFL) "\t" else "   ", row.names = FALSE, col.names = FALSE)
  
  write.table(c("", "> data import & export"),
              logCon, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  settings_table <- data.frame(
    col1 = c("nrows", "test_nrows", "header_translations", "comment.char", "column_separators", "selected separator", "na.strings", "header", "", "", ""),
    col3 = c(nrows, nrows_test, settings_header_input, encodeString(paste0("'", comment.char, "'")), encodeString(paste0("'", column_separators, "'", collapse = ", ")), encodeString(paste0("'", loadI$sep, "'")), paste(na.strings, collapse = ", "), header, "", "", ""),
    col4 = "   ",
    col5 = c("save_final_dataset", "order_columns", "out_header", "out_quote", "out_sep", "out_eol", "out_na", "out_colnames", "out_rownames", "out_dec", "out_qmethod"),
    col7 = c(save_final_dataset, order_columns, settings_header_output, out_quote, encodeString(paste0("'", out_sep, "'")), encodeString(paste0("'", out_eol, "'")), out_na, paste(out_colnames, collapse = ", "), paste(out_rownames, collapse = ", "), out_dec, out_qmethod),
    stringsAsFactors = FALSE)
  write.table(if(SFL) settings_table else format(settings_table, justify = "left"),
              logCon, quote = FALSE, sep = if(SFL) "\t" else "   ", row.names = FALSE, col.names = FALSE)
  close(logCon)
  rm(settings_table)
  

  # creating output file
  
  effect_quantile <- quantile(dataI$EFFECT, names = FALSE)
  
  useStrand <- strand_ref$plus + strand_ref$minus > 0L
  
  empty_colls <- c("CHR","POSITION","STRAND","PVALUE",
                   "EFF_ALL_FREQ","HWE_PVAL","CALLRATE","N_TOTAL","IMPUTED",
                   "USED_FOR_IMP", "IMP_QUALITY")[c(
                     if(iNA_chr_N + inv_chr_N == SNPn_preQC) TRUE else all(iNA$chr | inv$chr),
                     if(iNA_pos_N + inv_pos_N == SNPn_preQC) TRUE else all(iNA$pos | inv$pos),
                     strand_post$plus + strand_post$minus == 0L,
                     if(calculate_missing_p) FALSE else { if(iNA_p_N + inv_p_N == SNPn_preQC) TRUE else all(iNA$p | inv$p) },
                     if(useFRQ) FALSE else all(iNA$FRQ | inv$FRQ),
                     if(iNA_HWE_N + inv_HWE_N == SNPn_preQC) TRUE else all(iNA$HWE | inv$HWE),
                     if(iNA_cal_N + inv_cal_N == SNPn_preQC) TRUE else all(iNA$cal | inv$cal),
                     !useN,
                     SNPn_postQC_geno + SNPn_postQC_imp == 0L,
                     all(is.na(dataI$USED_FOR_IMP)),
                     if(iNA_impQ_N + inv_impQ_N == SNPn_preQC) TRUE else all(iNA$impQ | inv$impQ)
                   )	]
  
  
  QC_results <- list(QC_successful = TRUE,
                     filename_input = filename_input,
                     filename_output = paste0(filename, ".txt"),
                     
                     sample_size		= stat_N_max,
                     sample_size_HQ = stat_N_HQ,
                     lambda	= plot_output$lambda[1],
                     lambda_geno = plot_output$lambda[2],
                     lambda_imp	= plot_output$lambda[3],
                     
                     SNP_N_input			= SNPn_input,
                     SNP_N_input_monomorphic = monomorp_N,
                     SNP_N_input_monomorphic_identic_alleles = same_al_N,
                     SNP_N_input_chr		= if(remove_X | remove_Y | remove_XY | remove_M) sum(chr_X_N, chr_Y_N, chr_XY_N, chr_M_N, na.rm = TRUE) else NA,
                     
                     SNP_N_preQC			= SNPn_preQC,
                     SNP_N_preQC_unusable	= remove_L1_N,
                     SNP_N_preQC_invalid	= SNPn_preQC_inv,
                     SNP_N_preQC_min		= if(strand_pre$plus + strand_pre$minus == 0L) { NA } else { strand_pre$minus },
                     
                     SNP_N_midQC			= SNPn_ref,
                     SNP_N_midQC_min		= if(useStrand)				 { strand_ref$minus } else { NA },
                     SNP_N_midQC_min_std	= if(use_allele_std & useStrand) { if(strand_ref$minus == 0L) { 0L } else { allele_out_std$n_negative_strand } } else { NA },
                     SNP_N_midQC_min_alt	= if(use_allele_alt & useStrand) { if(strand_ref$minus == 0L) { 0L } else { allele_out_alt$n_negative_strand } } else { NA },
                     SNP_N_midQC_min_new	= if(				 useStrand) { if(strand_ref$minus == 0L) { 0L } else { allele_out_new$n_negative_strand } } else { NA },
                     
                     SNP_N_midQC_strandswitch_std		= if(use_allele_std)			 { allele_out_std$n_strandswitch } else { NA },
                     SNP_N_midQC_strandswitch_std_min	= if(use_allele_std & useStrand) { if(strand_ref$minus == 0L) { 0L } else { allele_out_std$n_negative_switch } } else { NA },
                     SNP_N_midQC_strandswitch_alt		= if(use_allele_alt)			 { allele_out_alt$n_strandswitch } else { NA },
                     SNP_N_midQC_strandswitch_alt_min	= if(use_allele_alt & useStrand) { if(strand_ref$minus == 0L) { 0L } else { allele_out_alt$n_negative_switch } } else { NA },
                     
                     SNP_N_midQC_mismatch		= N_mismatch,
                     SNP_N_midQC_mismatch_std	= if(use_allele_std)			{ allele_out_std$n_mismatch } else { NA },
                     SNP_N_midQC_mismatch_std_min	= if(use_allele_std & useStrand)	{ if(strand_ref$minus == 0L) { 0L } else { allele_out_std$n_negative_mismatch } } else { NA },
                     SNP_N_midQC_mismatch_alt	= if(use_allele_alt)			{ allele_out_alt$n_mismatch } else { NA },
                     SNP_N_midQC_mismatch_alt_min	= if(use_allele_alt & useStrand)	{ if(strand_ref$minus == 0L) { 0L } else { allele_out_alt$n_negative_mismatch } } else { NA },
                     
                     SNP_N_midQC_flip_std	= if(use_allele_std) { allele_out_std$n_flipped } else { NA },
                     SNP_N_midQC_flip_alt	= if(use_allele_alt) { allele_out_alt$n_flipped } else { NA },
                     SNP_N_midQC_flip_new	=							allele_out_new$n_flipped,
                     
                     SNP_N_midQC_ambiguous		= N_ambiguous,
                     SNP_N_midQC_ambiguous_std	= if(use_allele_std) { allele_out_std$n_ambiguous } else { NA },
                     SNP_N_midQC_ambiguous_alt	= if(use_allele_alt) { allele_out_alt$n_ambiguous } else { NA },
                     SNP_N_midQC_ambiguous_new	= allele_out_new$n_ambiguous,
                     
                     SNP_N_midQC_suspect	= N_suspect,
                     SNP_N_midQC_suspect_std	= if(use_allele_std) { allele_out_std$n_suspect } else { NA },
                     SNP_N_midQC_suspect_alt	= if(use_allele_alt) { allele_out_alt$n_suspect } else { NA },
                     
                     SNP_N_midQC_diffEAF	= N_diffEAF,
                     SNP_N_midQC_diffEAF_std	= if(use_allele_std) { allele_out_std$n_diffEAF } else { NA },
                     SNP_N_midQC_diffEAF_alt	= if(use_allele_alt) { allele_out_alt$n_diffEAF } else { NA },
                     
                     SNP_N_postQC		= SNPn_postQC,
                     SNP_N_postQC_geno		= SNPn_postQC_geno,
                     SNP_N_postQC_imp		= SNPn_postQC_imp,
                     SNP_N_postQC_invalid	= SNPn_postQC_inv,
                     SNP_N_postQC_min		= if(strand_post$plus + strand_post$minus == 0L) { NA } else { strand_post$minus },
                     SNP_N_postQC_HQ		= SNPn_postQC_HQ,
                     
                     fixed_HWE = if(ignore_impstatus) {
                       if(all(iNA$HWE | inv$HWE)) { "insuf. data"
                       } else { max(dataI$HWE_PVAL[!iNA$HWE & !inv$HWE]) == min(dataI$HWE_PVAL[!iNA$HWE & !inv$HWE]) }
                     } else {
                       if(all(iNA$HWE | inv$HWE | !geno_list)) { "insuf. data"
                       } else { max(dataI$HWE_PVAL[geno_list & !iNA$HWE & !inv$HWE]) == min(dataI$HWE_PVAL[geno_list & !iNA$HWE & !inv$HWE]) }
                     },
                     fixed_callrate = if(ignore_impstatus) {
                       if(all(iNA$cal | inv$cal)) { "insuf. data"
                       } else { max(dataI$CALLRATE[!iNA$cal & !inv$cal]) == min(dataI$CALLRATE[!iNA$cal & !inv$cal]) }
                     } else {
                       if(all(iNA$cal | inv$cal | !geno_list)) { "insuf. data"
                       } else { max(dataI$CALLRATE[geno_list & !iNA$cal & !inv$cal]) == min(dataI$CALLRATE[geno_list & !iNA$cal & !inv$cal]) }
                     },
                     fixed_sampleN  = if(useN) { stat_N_max == min(dataI$N_TOTAL, na.rm = iNA_N_N+inv_N_N>0L) } else "no data",
                     fixed_impQ = if(ignore_impstatus) {
                       if(all(iNA$impQ | inv$impQ)) { "insuf. data"
                       } else { max(dataI$IMP_QUALITY[!iNA$impQ & !inv$impQ]) == min(dataI$IMP_QUALITY[!iNA$impQ & !inv$impQ]) }
                     } else {
                       if(all(iNA$impQ | inv$impQ | !imp_list)) { "insuf. data"
                       } else { max(dataI$IMP_QUALITY[imp_list & !iNA$impQ & !inv$impQ]) == min(dataI$IMP_QUALITY[imp_list & !iNA$impQ & !inv$impQ]) }
                     },
                     
                     effect_25	= effect_quantile[2],
                     effect_mean	= mean(dataI$EFFECT),
                     effect_median = effect_quantile[3],
                     effect_75	= effect_quantile[4],
                     SE_median	= stat_SE,
                     SE_median_HQ= stat_SE_HQ,
                     skewness	= stat_skewness,
                     skewness_HQ = stat_skewness_HQ,
                     kurtosis	= stat_kurtosis,
                     kurtosis_HQ = stat_kurtosis_HQ,
                     
                     all_ref_std_name	= allele_name_std,
                     all_ref_alt_name	= allele_name_alt,
                     
                     all_MAF_std_r = if(use_allele_std) allele_out_std$FRQ_cor else NA,
                     all_MAF_alt_r = if(use_allele_alt) allele_out_alt$FRQ_cor else NA,
                     all_ambiguous_MAF_std_r = if(use_allele_std)	allele_out_std$FRQ_cor_ambiguous else NA,
                     all_ambiguous_MAF_alt_r = if(use_allele_alt)	allele_out_alt$FRQ_cor_ambiguous else NA,
                     all_non_ambig_MAF_std_r = if(use_allele_std) allele_out_std$FRQ_cor_nonambi else NA,
                     all_non_ambig_MAF_alt_r = if(use_allele_alt) allele_out_alt$FRQ_cor_nonambi else NA,
                     all_ref_changed = allele_ref_changed,
                     
                     effectsize_return = outcome_ES$return_ES,
                     effectsizes_HQ	= if(outcome_ES$return_ES) outcome_ES$HQ_effectsizes else NULL,
                     pvalue_r		= outcome_P,
                     visschers_stat	= stat_Visscher,
                     visschers_stat_HQ = stat_Vissc_HQ,
                     
                     columns_std_missing	= if(header_info$missing_N == 0L) 0L else paste(header_info$missing_h, collapse = ", "),
                     columns_std_empty		= if(length(empty_colls) == 0L)   0L else paste(empty_colls, collapse = ", "),
                     columns_unidentified	= if(header_info$unknown_N == 0L) 0L else paste(header_info$unknown_h, collapse = ", "),
                     
                     outcome_useFRQ	= useFRQ,
                     outcome_useHWE	= useHWE,
                     outcome_useCal	= useCal,
                     outcome_useImp	= useImp,
                     outcome_useMan	= useMan,
                     
                     settings_ignore_impstatus = ignore_impstatus,
                     settings_filter_NA_FRQ = NAfilter_FRQ,
                     settings_filter_NA_HWE = NAfilter_HWE,
                     settings_filter_NA_cal = NAfilter_cal,
                     settings_filter_NA_imp = NAfilter_imp,
                     settings_filter_HQ_FRQ = HQfilter_FRQ,
                     settings_filter_HQ_HWE = HQfilter_HWE,
                     settings_filter_HQ_cal = HQfilter_cal,
                     settings_filter_HQ_imp = HQfilter_imp,
                     settings_filter_QQ_FRQ = QQfilter_FRQ,
                     settings_filter_QQ_HWE = QQfilter_HWE,
                     settings_filter_QQ_cal = QQfilter_cal,
                     settings_filter_QQ_imp = QQfilter_imp
  )
  
  if(save_final_dataset) {
    print(" - saving cleaned dataset", quote = FALSE)
    flush.console()
    
    if(settings_header_output != "standard") {
      if(settings_header_output == "original") {
        if(order_columns) {
          header_new <- colnames(dataI)
          order_orig <- integer(length = header_info$header_N)
          for(hI in 1:header_info$header_N) { order_orig[hI] <- which(header_new == header_info$header_h[hI]) }
          colnames(dataI)[order_orig] <- header_orig
        } else {
          colnames(dataI)[1:length(header_orig)] <- header_orig
        }
      } else {
        header_out <- translate_header(header = colnames(dataI), standard = out_header[ ,1], alternative = out_header)
        if(header_out$missing_N > 0L) save_log(5, "saving file", "missing column", SNPn_postQC, SNPn_postQC, "-", paste("Unable to translate column(s)", paste(header_out$missing_h, collapse = ", ")), filename_dir)
        colnames(dataI) <- header_out$header_h
        
        if(settings_header_output == "GenABEL") {
          if(!"build" %in% header_out$unknown_h) dataI$build <- as.factor("unknown")
          dataI$effallele <- as.factor(dataI$allele1)
          if(!"pgc" %in% header_out$unknown_h) {
            dataI$pgc <- if(is.na(plot_output$lambda[1])) NA else {
              pchisq( (dataI$beta/(dataI$sebeta * sqrt(plot_output$lambda[1]) ) )^2,
                      1, lower.tail=FALSE) }
          }
          if(!"lambda.estimate" %in% header_out$unknown_h) dataI$lambda.estimate <- plot_output$lambda[1]
          if(!"lambda.se" %in% header_out$unknown_h) dataI$lambda.se <- NA
          
          dataI <- dataI[ , c("name", "chromosome", "position", "strand",
                              "allele1", "allele2", "build", "effallele",
                              "effallelefreq", "n", "beta", "sebeta", "p",
                              "pgc", "lambda.estimate", "lambda.se",
                              "pexhwe", "call")]
        }
      }
    }
    
    write.table(dataI,
                if(gzip_final_dataset) gzfile(paste0(filename_dir, ".txt.gz")) else paste0(filename_dir, ".txt"),
                quote = out_quote, sep = out_sep, eol = out_eol, na = out_na, dec = out_dec,
                row.names = out_rownames, col.names = out_colnames, qmethod = out_qmethod)
  }
  
  print("", quote = FALSE)
  print(paste("QC check completed for", filename_input,"( file", logI, "out of", logN,")"), quote = FALSE)
  return(QC_results)
}
