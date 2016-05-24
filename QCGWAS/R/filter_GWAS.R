filter_GWAS <-
function(ini_file,
                        GWAS_files,
                        output_names,
                        gzip_output = TRUE,
                        
                        dir_GWAS = getwd(),
                        dir_output = dir_GWAS,
                        
                        FRQ_HQ = NULL, HWE_HQ = NULL, cal_HQ = NULL, imp_HQ = NULL,
                        FRQ_NA = TRUE, HWE_NA = TRUE, cal_NA = TRUE, imp_NA = TRUE,
                        ignore_impstatus = FALSE,
                        
                        remove_X = FALSE, remove_Y = FALSE, remove_XY = FALSE, remove_M = FALSE,
                        
                        header_translations,
                        check_impstatus = FALSE,
                        imputed_T = c("1", "TRUE", "yes", "YES", "y", "Y"), imputed_F = c("0", "FALSE", "no", "NO", "n", "N"), imputed_NA = NULL,
                        
                        column_separators = c("\t", " ", "", ",", ";"), header = TRUE,
                        nrows = -1, nrows_test = 1000, comment.char = "",
                        na.strings = c("NA", "."),
                        
                        out_header = "original", out_quote = FALSE, out_sep = "\t",
                        out_eol = "\n", out_na = "NA", out_dec = ".", out_qmethod = "escape",
                        out_rownames = FALSE, out_colnames = TRUE, ...){
  
  # Checking input arguments
  
  stopifnot(is.character(dir_GWAS), length(dir_GWAS) == 1L, is.character(dir_output), length(dir_output) == 1L)
  if(!file.exists(dir_GWAS)) stop("Cannot find folder 'dir_GWAS'")
  if(!file.exists(dir_output)) stop("Cannot find folder 'dir_output'")
  stopifnot(is.logical(gzip_output), length(gzip_output) == 1L)
  if(is.na(gzip_output)) stop("'gzip_output' cannot be NA")
  stopifnot(is.logical(header), length(header) == 1L)
  if(is.na(header)) stop("'header' cannot be NA")
  stopifnot(is.numeric(nrows), length(nrows) == 1L)
  if(is.na(nrows) | nrows == 0) stop("'nrows' cannot be missing or zero")
  stopifnot(is.numeric(nrows_test), length(nrows_test) == 1L)
  if(is.na(nrows_test) | nrows_test == 0) stop("'nrows_test' cannot be missing or zero")
  stopifnot(is.character(comment.char), length(comment.char) == 1L)
  if(nchar(comment.char) != 1L & nchar(comment.char) != 0L) stop("'comment.char' must be a single character or empty string")
  
  stopifnot(is.character(na.strings), is.character(column_separators),
            is.logical(check_impstatus), length(check_impstatus) == 1L,
            !is.na(check_impstatus))
  if(check_impstatus){
    stopifnot(is.character(imputed_T), is.character(imputed_F), is.character(imputed_NA))
    if(any(duplicated(c(imputed_T, imputed_F, imputed_NA)))) stop("duplicate strings in the 'imputed' arguments")
  }
  
  stopifnot(is.logical(out_quote), length(out_quote) == 1L,
            is.character(out_sep), length(out_sep) == 1L,
            is.character(out_eol), length(out_eol) == 1L,
            is.character(out_na),  length(out_na) == 1L,
            is.character(out_dec), length(out_dec) == 1L,
            is.character(out_qmethod), length(out_qmethod) == 1L)
  if(is.na(out_quote)) stop("'out_quote' cannot be NA")
  if(nchar(out_sep) == 0L) stop("'out_sep' cannot be an empty character-string")
  if(nchar(out_eol) == 0L) stop("'out_eol' cannot be an empty character-string")
  if(nchar(out_dec) != 1L) stop("'out_dec' must be of length 1")
  if(!out_qmethod %in% c("escape", "double", "e", "d")) stop("'out_qmethod' must be either 'escape' or 'double'")
  
  stopifnot(is.logical(remove_X), length(remove_X) == 1L,
            is.logical(remove_Y), length(remove_Y) == 1L,
            is.logical(remove_XY), length(remove_XY) == 1L,
            is.logical(remove_M), length(remove_M) == 1L)
  if(is.na(remove_X)) stop( "'remove_X' cannot be NA")
  if(is.na(remove_Y)) stop( "'remove_Y' cannot be NA")
  if(is.na(remove_XY))stop("'remove_XY' cannot be NA")
  if(is.na(remove_M)) stop( "'remove_M' cannot be NA")
  
  if(missing(header_translations)) {    
    check_header <- FALSE  
  } else {
    if(is.character(header_translations) & !is.matrix(header_translations)) {
      if(length(header_translations) != 1L & nchar(header_translations) < 3L) {
        stop("'header_translations' is neither a table nor a filename") }
      print(paste("Loading 'header_translations' from file:", header_translations), quote = FALSE)
      if(file.exists(paste(dir_GWAS, header_translations, sep = "/"))) {
        header_translations <- read.table(paste(dir_GWAS, header_translations, sep = "/"), stringsAsFactors = FALSE)
      } else {
        if(file.exists(header_translations)) {
          header_translations <- read.table(header_translations, stringsAsFactors = FALSE)
        } else { stop("Cannot find 'header_translations'")}
      }
    }
    if(!is.data.frame(header_translations) & !is.matrix(header_translations)) stop("'header_translations' is not a table")
    if(ncol(header_translations) != 2L) stop("'header_translations' does not have two columns")
    if(any(duplicated(header_translations[,2]))) stop("'header_translations' contains duplicated elements in column 2")
    check_header <- TRUE
  }
  
  out_header_type <- "table"
  if(is.character(out_header) & !is.matrix(out_header)) {
    if(length(out_header) != 1L | nchar(out_header) < 3L) stop("'out_header' is neither a table nor a filename")
    if(out_header %in% c("original", "standard", "GWAMA", "PLINK", "META", "GenABEL", "old")) {
      out_header_type <- out_header
      if(out_header_type != "standard" & out_header_type != "original") {
        if(out_header_type == "GenABEL"){
          out_header <- data.frame(
            GenABEL = c("name", "chromosome", "position", "strand", "allele1", "allele2", "effallelefreq", "n", "beta", "sebeta", "p", "pexhwe", "call"),
            QC = c("MARKER", "CHR", "POSITION", "STRAND", "EFFECT_ALL", "OTHER_ALL", "EFF_ALL_FREQ", "N_TOTAL", "EFFECT", "STDERR", "PVALUE", "HWE_PVAL", "CALLRATE"),
            stringsAsFactors = FALSE)
        } else {
          out_header <- data.frame(
            GWAMA = c("MARKER", "CHR", "POSITION", "EA", "NEA", "STRAND", "BETA", "SE", "P", "EAF", "N", "IMPUTED", "IMP_QUALITY"),
            PLINK = c("SNP",    "CHR", "BP",  		 "A1", "A2",	"STRAND", "BETA", "SE", "P", "EFF_ALL_FREQ", "N", "IMPUTED", "IMP_QUALITY"),
            META = c( "rsid",	 "chr", "pos",			"allele_B", "allele_A", "strand", "beta", "se", "P_value", "EFF_ALL_FREQ", "N", "imputed", "info"),
            old =c("MARKER", "CHR", "POSITION", "ALLELE1",		"ALLELE2",	 "STRAND", "EFFECT", "STDERR", "PVALUE", "FREQLABEL",		"N_TOTAL", "IMPUTED", "IMP_QUALITY"),
            QC = c("MARKER", "CHR", "POSITION", "EFFECT_ALL", "OTHER_ALL", "STRAND", "EFFECT", "STDERR", "PVALUE", "EFF_ALL_FREQ", "N_TOTAL", "IMPUTED", "IMP_QUALITY"),
            stringsAsFactors = FALSE)
          out_header <- out_header[ ,c(out_header_type, "QC")]
        }
      }
    } else {
      print(paste("Loading 'out_header' from file:", out_header), quote = FALSE)
      if(file.exists(paste(dir_GWAS, out_header, sep = "/"))) {
        out_header <- read.table(paste(dir_GWAS, out_header, sep = "/"), stringsAsFactors = FALSE)
      } else {
        if(file.exists(out_header)) {
          out_header <- read.table(out_header, stringsAsFactors = FALSE)
        } else {stop("'out_header' is neither a standard type nor a file name")}}
    }  
  }
  if(out_header_type == "table") {
    if(is.matrix(out_header) | is.data.frame(out_header)) {
      if(ncol(out_header) != 2L) stop("'out_header' does not have 2 columns")
      if(any(is.na(out_header))) stop("'out_header' contains missing values")
      if(any(duplicated(out_header[,1]),duplicated(out_header[,2]))) stop("'out_header' contains duplicated names")
    } else { stop("'out_header' is not a table or a data-frame") } }
  
  
  # Checking the ini file. If not present, a table is 
  # created from the other input arguments (tested here)
  
  if(missing(ini_file)){
    if(missing(GWAS_files)) {
      ini_check <- TRUE
      if(file.exists(paste(dir_GWAS, "Check_filtersettings.txt", sep = "/"))) {
        print("No ini-file specified - loading 'Check_filtersettings.txt'", quote = FALSE)
        ini_file <- read.table(paste(dir_GWAS, "Check_filtersettings.txt", sep = "/"),
                               header = TRUE, stringsAsFactors = FALSE)
      } else {
        if(file.exists("Check_filtersettings.txt")) {
          print("No ini-file specified - loading 'Check_filtersettings.txt'", quote = FALSE)
          ini_file <- read.table("Check_filtersettings.txt",
                                 header = TRUE, stringsAsFactors = FALSE)
        } else { stop("No ini_file specified") }
      }
    } else {
      ini_check <- FALSE
      if(!is.character(GWAS_files)) stop("'GWAS_files' is not a character string or vector")
      if(length(GWAS_files) == 1L) GWAS_files <- list.files(path = dir_GWAS, pattern = GWAS_files)
      file_N <- length(GWAS_files)
      if(file_N == 0L) stop("'GWAS_files' does not point to existing files")
      
      stopifnot(is.logical(FRQ_NA),
                is.logical(HWE_NA),
                is.logical(cal_NA),
                is.logical(imp_NA),
                is.logical(ignore_impstatus))
      if(is.matrix(FRQ_NA) | is.matrix(HWE_NA) | is.matrix(cal_NA) | is.matrix(imp_NA) | is.na(ignore_impstatus)) {
        stop("'xxx_NA' and 'ignore_impstatus' arguments cannot be matrices") }
      
      if(length(FRQ_NA) != 1L & length(FRQ_NA) != file_N) {
        if(file_N %% length(FRQ_NA) == 0L & length(FRQ_NA) < file_N) {
          print("Warning: vector 'FRQ_NA' is recycled to bring it to the same length as the file list")
        } else { stop("The length of 'FRQ_NA' does not match the file-list") } }
      if(length(HWE_NA) != 1L & length(HWE_NA) != file_N) {
        if(file_N %% length(HWE_NA) == 0L & length(HWE_NA) < file_N) {
          print("Warning: vector 'HWE_NA' is recycled to bring it to the same length as the file list")
        } else { stop("The length of 'HWE_NA' does not match the file-list") } }
      if(length(cal_NA) != 1L & length(cal_NA) != file_N) {
        if(file_N %% length(cal_NA) == 0L & length(cal_NA) < file_N) {
          print("Warning: vector 'cal_NA' is recycled to bring it to the same length as the file list")
        } else { stop("The length of 'cal_NA' does not match the file-list") } }
      if(length(imp_NA) != 1L & length(imp_NA) != file_N) {
        if(file_N %% length(imp_NA) == 0L & length(imp_NA) < file_N) {
          print("Warning: vector 'imp_NA' is recycled to bring it to the same length as the file list")
        } else { stop("The length of 'imp_NA' does not match the file-list") } }
      if(length(ignore_impstatus) != 1L & length(ignore_impstatus) != file_N) {
        if(file_N %% length(ignore_impstatus) == 0L & length(ignore_impstatus) < file_N) {
          print("Warning: vector 'ignore_impstatus' is recycled to bring it to the same length as the file list")
        } else { stop("The length of 'ignore_impstatus' does not match the file-list") } }
        
      if(is.null(FRQ_HQ)) {
        FRQ_HQ <- NA
        FRQ_NA <- FALSE
      } else {
        if(!is.vector(FRQ_HQ)) stop("'FRQ_HQ' must be a vector")
        if(!is.numeric(FRQ_HQ) & !all(is.na(FRQ_HQ))) stop("'FRQ_HQ' isn't a numerical or NA value") }
      if(is.null(HWE_HQ)) {
        HWE_HQ <- NA
        HWE_NA <- FALSE
      } else {
        if(!is.vector(HWE_HQ)) stop("'HWE_HQ' must be a vector")
        if(!is.numeric(HWE_HQ) & !all(is.na(HWE_HQ))) stop("'HWE_HQ' isn't a numerical or NA value") }
      if(is.null(cal_HQ)) {
        cal_HQ <- NA
        cal_NA <- FALSE
      } else {
        if(!is.vector(cal_HQ)) stop("'cal_HQ' must be a vector")
        if(!is.numeric(cal_HQ) & !all(is.na(cal_HQ))) stop("'cal_HQ' isn't a numerical or NA value") }
      if(is.null(imp_HQ)) {
        imp_HQ <- NA
        imp_NA <- FALSE
      } else {
        if(!is.vector(imp_HQ)) stop("'imp_HQ' must be a vector")
        if(!is.numeric(imp_HQ) & !all(is.na(imp_HQ))) stop("'imp_HQ' isn't a numerical or NA value") }
      
      if(length(FRQ_HQ) != 1L & length(FRQ_HQ) != file_N) {
        if(file_N %% length(FRQ_HQ) == 0L & length(FRQ_HQ) < file_N) {
          print("Warning: vector 'FRQ_HQ' is recycled to bring it to the same length as the file list")
        } else { stop("The length of 'FRQ_HQ' does not match the file-list") } }
      if(length(HWE_HQ) != 1L & length(HWE_HQ) != file_N) {
        if(file_N %% length(HWE_HQ) == 0L & length(HWE_HQ) < file_N) {
          print("Warning: vector 'HWE_HQ' is recycled to bring it to the same length as the file list")
        } else { stop("The length of 'HWE_HQ' does not match the file-list") } }
      if(length(cal_HQ) != 1L & length(cal_HQ) != file_N) {
        if(file_N %% length(cal_HQ) == 0L & length(cal_HQ) < file_N) {
          print("Warning: vector 'cal_HQ' is recycled to bring it to the same length as the file list")
        } else { stop("The length of 'cal_HQ' does not match the file-list") } }
      if(length(imp_HQ) != 1L & length(imp_HQ) != file_N) {
        if(file_N %% length(imp_HQ) == 0L & length(imp_HQ) < file_N) {
          print("Warning: vector 'imp_HQ' is recycled to bring it to the same length as the file list")
        } else { stop("The length of 'imp_HQ' does not match the file-list") } }
      
      ini_file <- data.frame(file = GWAS_files,
                             filter_FRQ = FRQ_HQ, filter_HWE = HWE_HQ,
                             filter_cal = cal_HQ, filter_imp = imp_HQ,
                             FRQ_NA = FRQ_NA, HWE_NA = HWE_NA,
                             cal_NA = cal_NA, imp_NA = imp_NA,
                             ignore_impstatus = ignore_impstatus,
                             out_name = "",
                             stringsAsFactors = FALSE)
      
      if(any(is.na(ini_file$FRQ_NA)) | any(is.na(ini_file$HWE_NA)) |
           any(is.na(ini_file$cal_NA)) | any(is.na(ini_file$imp_NA)) |
           any(is.na(ini_file$ignore_impstatus))) {
        stop("'xxx_NA' and 'ignore_impstatus' arguments cannot be NA") } }
  } else {
    ini_check <- TRUE
    if(is.character(ini_file)){
      if(file.exists(paste(dir_GWAS, ini_file, sep = "/"))) {
        ini_file <- read.table(paste(dir_GWAS, ini_file, sep = "/"),
                               header = TRUE, stringsAsFactors = FALSE)
      } else {
        if(file.exists(ini_file)) {
          ini_file <- read.table(ini_file,
                                 header = TRUE, stringsAsFactors = FALSE)
        } else { stop("Cannot find the specified ini file") }
  } } }
  
  # Checking the contents of the ini file (if present)
  if(ini_check) {
    stopifnot(is.data.frame(ini_file), nrow(ini_file) > 0L,
              all(c("file",
                    "filter_FRQ", "filter_HWE",
                    "filter_cal", "filter_imp",
                    "FRQ_NA", "HWE_NA",
                    "cal_NA", "imp_NA", "ignore_impstatus") %in% colnames(ini_file)))
    
    stopifnot(is.numeric(ini_file$filter_FRQ) | all(is.na(ini_file$filter_FRQ)),
              is.numeric(ini_file$filter_HWE) | all(is.na(ini_file$filter_HWE)),
              is.numeric(ini_file$filter_cal) | all(is.na(ini_file$filter_cal)),
              is.numeric(ini_file$filter_imp) | all(is.na(ini_file$filter_imp)),
              is.logical(ini_file$FRQ_NA) & !any(is.na(ini_file$FRQ_NA)),
              is.logical(ini_file$HWE_NA) & !any(is.na(ini_file$HWE_NA)),
              is.logical(ini_file$cal_NA) & !any(is.na(ini_file$cal_NA)),
              is.logical(ini_file$imp_NA) & !any(is.na(ini_file$imp_NA)),
              is.logical(ini_file$ignore_impstatus) & !any(is.na(ini_file$ignore_impstatus)))
    file_N <- nrow(ini_file)
  }
  
  if(!all(file.exists(paste(dir_GWAS, ini_file$file, sep = "/")))) {
    stop(paste("Cannot find file(s):",
               paste(ini_file$file[!file.exists(paste(dir_GWAS, ini_file$file, sep = "/"))],
                     collapse = ", "))) }
  
  if(missing(output_names)){
    ini_file$out_name <- ini_file$file
  } else {
    stopifnot(is.character(output_names), length(output_names) == file_N, all(nchar(output_names) > 2L))
    ini_file$out_name <- output_names
  }
  
  
  # The actual filtering starts here
  print("", quote = FALSE)
  print("Processing:", quote = FALSE)
  print(paste(" -", 1:file_N, ini_file$file), quote = FALSE)
  print("", quote = FALSE)
  
  return_list <- logical(length = file_N)
  use_Chr <- remove_X | remove_Y | remove_XY | remove_M
  
  for(fi in 1:file_N){
    print("", quote = FALSE)
    print(paste0("Loading ", ini_file$file[fi], " - (", fi, "/", file_N, ")"), quote = FALSE)
    flush.console()
    dataI <- load_GWAS(filename = ini_file$file[fi], dir = dir_GWAS,
                       column_separators = column_separators, test_nrows = nrows_test,
                       header = header, nrows = nrows,  comment.char = comment.char,
                       na.strings = na.strings, ...)
    
    print(paste("Processing", ini_file$file[fi]), quote = FALSE)
    flush.console()
    
    use_FRQ <- ini_file$FRQ_NA[fi] | !is.na(ini_file$filter_FRQ[fi])
    use_HWE <- ini_file$HWE_NA[fi] | !is.na(ini_file$filter_HWE[fi])
    use_cal <- ini_file$cal_NA[fi] | !is.na(ini_file$filter_cal[fi])
    use_imp <- ini_file$imp_NA[fi] | !is.na(ini_file$filter_imp[fi])
    
    use_impstatus <- if(ini_file$ignore_impstatus[fi]) check_impstatus else use_HWE | use_cal | use_imp | check_impstatus
    required_cols <- c("CHR","EFF_ALL_FREQ","HWE_PVAL","CALLRATE", "IMP_QUALITY", "IMPUTED")[
      c(use_Chr, use_FRQ, use_HWE, use_cal, use_imp, use_impstatus)]
    
    # Checking the file header
    if(check_header){
      
      header_orig <- colnames(dataI)
      header_info <- translate_header(header = header_orig, alternative = header_translations)
      
      if(any(duplicated(header_info$header_h))) { stop(paste("Dataset contains duplicate columns:", paste(header_info$header_h[duplicated(header_info$header_h)], collapse = ", "))) }
      if(header_info$missing_N > 0L) {
        if(any(required_cols %in% header_info$missing_h)) {
          print(paste(" - ERROR: column(s)", paste(required_cols[required_cols %in% header_info$missing_h], collapse = ", "), "not found: cannot apply filters"), quote = FALSE)
          next
      } }
      colnames(dataI) <- header_info$header_h
    } else {
      if(any(!required_cols %in% colnames(dataI))) {
        print(paste(" - ERROR: column(s)", paste(required_cols[!required_cols %in% colnames(dataI)], collapse = ", "), "not found: cannot apply filters"), quote = FALSE)
        next
    } }
    
    # Convert (if required) and check imputatution status
    if(use_impstatus){
      if(check_impstatus) { dataI$IMPUTED <- convert_impstatus(dataI$IMPUTED, T_strings = imputed_T, F_strings = imputed_F, NA_strings = imputed_NA, use_log = FALSE) }
      if(any(is.na(dataI$IMPUTED))){
        if(all(is.na(dataI$IMPUTED))) {
          print(" - ERROR: imputation status missing or untranslated. Dataset skipped.", quote = FALSE)
          next
        }
        print(" - WARNING: missing or untranslated imputation-status values", quote = FALSE)
    } }
    
    # Applying filter
    filter_list	<- HQ_filter(data = dataI, ignore_impstatus = ini_file$ignore_impstatus[fi],
                             FRQ_val = if(use_FRQ) ini_file$filter_FRQ[fi] else NULL,
                             HWE_val = if(use_HWE) ini_file$filter_HWE[fi] else NULL,
                             cal_val = if(use_cal) ini_file$filter_cal[fi] else NULL,
                             imp_val = if(use_imp) ini_file$filter_imp[fi] else NULL,
                             FRQ_NA  = ini_file$FRQ_NA[fi], HWE_NA = ini_file$HWE_NA[fi],
                             cal_NA  = ini_file$cal_NA[fi], imp_NA = ini_file$imp_NA[fi]) 
    
    if(use_Chr){
      notNA <- !is.na(dataI$CHR)
      na_rm_chr <- any(!notNA)
      if(remove_X) {
        dataI$CHR[dataI$CHR %in% c("X", "x")] <- 23L
        chr_N <- sum(dataI$CHR == 23, na.rm = na_rm_chr)
        print(paste(" - Chr. X SNPs removed:", chr_N), quote = FALSE)
        if(chr_N > 0L) { filter_list[dataI$CHR == 23 & notNA] <- FALSE }
      }
      if(remove_Y) {
        dataI$CHR[dataI$CHR %in% c("Y", "y")] <- 24L
        chr_N <- sum(dataI$CHR == 24, na.rm = na_rm_chr)
        print(paste(" - Chr. Y SNPs removed:", chr_N), quote = FALSE)
        if(chr_N > 0L) { filter_list[dataI$CHR == 24 & notNA] <- FALSE }
      }
      if(remove_XY) {
        dataI$CHR[dataI$CHR %in% c("XY","xy", "Xy", "xY")]<- 25L
        chr_N <- sum(dataI$CHR == 25, na.rm = na_rm_chr)
        print(paste(" - Chr. XY SNPs removed:", chr_N), quote = FALSE)
        if(chr_N > 0L) { filter_list[dataI$CHR == 25 & notNA] <- FALSE }
      }
      if(remove_M) {
        dataI$CHR[dataI$CHR %in% c("M", "m", "MT", "mt", "Mt", "mT")] <- 26L
        chr_N <- sum(dataI$CHR == 26, na.rm = na_rm_chr)
        print(paste(" - Chr. M SNPs removed:", chr_N), quote = FALSE)
        if(chr_N > 0L) { filter_list[dataI$CHR == 26 & notNA] <- FALSE }
      }
    }
    
    filter_N <- sum(!filter_list)
    if(filter_N == nrow(dataI)) {
      print(" - WARNING: no SNPs left after filtering", quote = FALSE)
    } else {
      if(filter_N > 0L) { dataI <- dataI[filter_list, ] } 
      print(paste(" - SNPs filtered:", filter_N), quote = FALSE)
      flush.console()
      
      # Changing header & saving data
      if(out_header_type != "standard") {
        if(out_header_type == "original") { if(check_header) { colnames(dataI) <- header_orig }
        } else {
          header_out <- translate_header(header = colnames(dataI),
                                         standard = out_header[,1], alternative = out_header)
          if(header_out$missing_N > 0L) { print(paste("Unable to translate column(s)", paste(header_out$missing_h, collapse = ", ")), quote = FALSE) }
          colnames(dataI) <- header_out$header_h

          if(out_header_type == "GenABEL") {
            if(!"build" %in% header_out$unknown_h) dataI$build <- as.factor("unknown")
            dataI$effallele <- as.factor(dataI$allele1)
            if(!"pgc" %in% header_out$unknown_h) dataI$pgc <- NA
            if(!"lambda.estimate" %in% header_out$unknown_h) dataI$lambda.estimate <- NA
            if(!"lambda.se" %in% header_out$unknown_h) dataI$lambda.se <- NA
            
            dataI <- dataI[ , c("name", "chromosome", "position", "strand",
                                "allele1", "allele2", "build", "effallele",
                                "effallelefreq", "n", "beta", "sebeta", "p",
                                "pgc", "lambda.estimate", "lambda.se",
                                "pexhwe", "call")]
          }
      } }
      
      write.table(dataI,
                  if(gzip_output) gzfile(paste0(dir_output, "/", ini_file$out_name[fi], ".gz")) else paste(dir_output, ini_file$out_name[fi], sep = "/"),
                  quote = out_quote, sep = out_sep, eol = out_eol, na = out_na, dec = out_dec,
                  row.names = out_rownames, col.names = out_colnames, qmethod = out_qmethod)
      return_list[fi] <- TRUE
    }
  } # end of for loop
  if(check_header) rm(header_translations)
  if(use_Chr) rm(notNA)
  rm(ini_file, out_header, dataI, filter_list)
  gc(FALSE)
  return(invisible(return_list))
}
