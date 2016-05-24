QC_plots <-
function(dataset, plot_QQ = TRUE, plot_Man = TRUE,
                     FRQfilter_values = NULL, FRQfilter_NA = filter_NA,
                     HWEfilter_values = NULL, HWEfilter_NA = filter_NA,
                     calfilter_values = NULL, calfilter_NA = filter_NA,
                     impfilter_values = NULL, impfilter_NA = filter_NA, impfilter_min = min(dataset$IMP_QUALITY, na.rm = TRUE),
                     manfilter_FRQ = NULL, manfilter_HWE = NULL, manfilter_cal = NULL, manfilter_imp = NULL,
                     filter_NA = TRUE, plot_cutoff_p = 0.05, plot_names = FALSE,
                     QQ_colors = c("red", "blue", "orange", "green3", "yellow"), plot_QQ_bands = FALSE,
                     save_name = "dataset", save_dir = getwd(), header_translations, use_log = FALSE,
                     check_impstatus = FALSE, ignore_impstatus = FALSE, T_strings = c("1", "TRUE", "yes", "YES", "y", "Y"), F_strings = c("0", "FALSE", "no", "NO", "n", "N"), NA_strings = c(NA, "NA", ".", "-")
) {
  
  if(!is.null(manfilter_FRQ)) { if(is.na(manfilter_FRQ) & !FRQfilter_NA) manfilter_FRQ <- NULL }
  if(!is.null(manfilter_HWE)) { if(is.na(manfilter_HWE) & !HWEfilter_NA) manfilter_HWE <- NULL }
  if(!is.null(manfilter_cal)) { if(is.na(manfilter_cal) & !calfilter_NA) manfilter_cal <- NULL }
  if(!is.null(manfilter_imp)) { if(is.na(manfilter_imp) & !impfilter_NA) manfilter_imp <- NULL }
  
  col_FRQ <- !is.null(FRQfilter_values) | (plot_Man & !is.null(manfilter_FRQ))
  col_HWE <- !is.null(HWEfilter_values) | (plot_Man & !is.null(manfilter_HWE))
  col_cal <- !is.null(calfilter_values) | (plot_Man & !is.null(manfilter_cal))
  col_impQ<- !is.null(impfilter_values) | (plot_Man & !is.null(manfilter_imp))
  col_N	  <- any(c(FRQfilter_values, HWEfilter_values, calfilter_values, impfilter_values) >= 1, na.rm = TRUE)
  col_impS<- TRUE
  
  if(is.vector(dataset)) { dataset <- data.frame(PVALUE = dataset) }
  dataN <- nrow(dataset)
  header_std <- c("MARKER", "PVALUE", "CHR", "POSITION",
                  "EFF_ALL_FREQ", "HWE_PVAL", "CALLRATE", "IMP_QUALITY",
                  "IMPUTED", "N_TOTAL")[c(plot_names, TRUE, plot_Man, plot_Man,
                                          col_FRQ, col_HWE, col_cal, col_impQ,
                                          col_impS, col_N)]
  
  check_header <- !missing(header_translations)
  if(check_header) {
    header_dat <- toupper(colnames(dataset))
    header_col <- integer(length = length(header_std))
    for(forI in 1:length(header_std)) {
      header_cur <- identify_column(header_std[forI], header_translations, header_dat)
      if(length(header_cur) == 1L) { header_col[forI] <- header_cur
      } else {
        if(length(header_cur) == 0L) {
          if(header_std[forI] == "IMPUTED") { col_impS <- FALSE
          } else { stop(paste("Cannot identify column:", header_std[forI])) }
        } else { stop(paste("Multiple columns identified as:", header_std[forI])) }
      }
    }
    if(!col_impS) { 
      header_std <- c("MARKER", "PVALUE", "CHR", "POSITION", "EFF_ALL_FREQ", "HWE_PVAL", "CALLRATE", "IMP_QUALITY", "IMPUTED", "N_TOTAL")[c(plot_names, TRUE, plot_Man, plot_Man, col_FRQ, col_HWE, col_cal, col_impQ, col_impS, col_N)]
      header_col <- header_col[!header_col == 0L]
    }
    dataset <- dataset[ , header_col]
    colnames(dataset) <- header_std
  } else {
    header_missing <- header_std[!header_std %in% colnames(dataset)]
    if(length(header_missing) > 0L) {
      if(length(header_missing) == 1L & header_missing[1] == "IMPUTED") { col_impS <- FALSE
      } else { stop(paste("Missing columns:", paste(header_missing, collapse = ", "))) }
    }
  }
  
  if(col_impS) {
    if(check_impstatus) {
      if(!check_header) dataset <- dataset[ , header_std]
      dataset$IMPUTED <- convert_impstatus(dataset$IMPUTED, T_strings = T_strings, F_strings = F_strings, NA_strings = NA_strings, use_log = FALSE)
    }
    no_impstatus <- all(is.na(dataset$IMPUTED))
    if(no_impstatus) {
      if (use_log) { save_log(phaseL = 4L, checkL = "QC statistics", typeL = "imputation status", SNPL = dataN, allSNPs = dataN, actionL = if(ignore_impstatus) "-" else "did not filter", noteL = "Cannot parse imputation-status column", fileL = paste(save_dir, save_name, sep = "/"))
      } else { print("No valid, non-missing imputation-status values", quote = FALSE) }
    }
  } else {
    if(check_impstatus) stop("Cannot check impstatus - no such column in dataset")
    print(" - - cannot not identify imputation-status column", quote = FALSE)
    no_impstatus <- TRUE
  }
  
  if(col_HWE | col_cal | col_impQ) {
    if(no_impstatus) {
      if(!ignore_impstatus) {
        if(!col_impQ) {
          print(" - - No imputation status: all SNPs set to genotyped", quote = FALSE)
          ignore_impstatus <- TRUE
        } else {
          if(!col_HWE & !col_cal) {
            print(" - - No imputation status: all SNPs set to imputed", quote = FALSE)
            ignore_impstatus <- TRUE
          } else stop(" - - No imputation status: cannot apply filters")
        }
      }
    } else {
      genoset <- dataset$IMPUTED == 0L & !is.na(dataset$IMPUTED)
      imp_set <- dataset$IMPUTED == 1L & !is.na(dataset$IMPUTED)
      
      if(!ignore_impstatus) {
        if((col_HWE | col_cal) & !any(genoset)) {
          print(" - - No genotyped SNPs: HWE p-value & call rate filters disabled", quote = FALSE)
          HWEfilter_values	<- NULL
          calfilter_values	<- NULL
          manfilter_HWE	<- NULL
          manfilter_cal	<- NULL
        }
        if(col_impQ & !any(imp_set)) {
          print(" - - No imputed SNPs found: imputation-quality filter disabled", quote = FALSE)
          impfilter_values	<- NULL
          manfilter_imp	<- NULL
      }	} }
  } else {
    if(!no_impstatus) {
      genoset <- dataset$IMPUTED == 0 & !is.na(dataset$IMPUTED)
      imp_set <- dataset$IMPUTED == 1 & !is.na(dataset$IMPUTED)
  } }
  
  
  # Stage 1: creating filters & calculating lambda
  if(!is.null(c(FRQfilter_values, HWEfilter_values, calfilter_values, impfilter_values))) {
    if(col_N) {	if(any(is.na(dataset$N_TOTAL))) {
      if (use_log) { save_log(phaseL = 4L, checkL = "QC statistics", typeL = "sample size", SNPL = sum(is.na(dataset$N_TOTAL)), allSNPs = dataN, actionL = "did not filter", noteL = "Missing sample-size: size-based filters will ignore or exclude these entries according to their filter_NA setting", fileL = paste(save_dir, save_name, sep = "/"))
      } else { print(paste(" - - Warning: missing sample sizes (", round(sum(is.na(dataset$N_TOTAL))/dataN, digits = 0), "% of SNPs) in dataset. Sample size filters will ignore or exclude these entries according to their respective filter_NA setting."), quote = FALSE) }
    } }
    flush.console()
    
    # creating functions for the QQ filter
    QQfilter_sort <- function(..., minimal_value = 0) {
      if(is.null(c(...))) stop("no QQ-filter values specified")
      if(!is.numeric(c(...))) stop("non-numeric QQ-filter input")
      if(minimal_value >= 1 ) stop("invalid minimal or maximal value")
      
      input <- unique(c(...))
      if(length(input) > 5L) input <- input[1:5]
      output <- numeric(length = 0)
      if(any(is.na(input) | input <= minimal_value)) {
        output <- c(output, minimal_value)
        if(any(is.na(input))) input[is.na(input)] <- minimal_value # removing NA's from input
      }
      if(any(input >= 1)) output <- c(output, sort(input[which(input >= 1)]) )
      if(any(input <  1 & input > minimal_value)) output <- c(output, sort(input[which(input < 1 & input > minimal_value)]) )
      return(output)
    }
    
    QQfilter_name <- function(var_name, filter_values, filter_NA = TRUE, minimal_value = 0, sign = ">=") {
      name_list <- character(length = length(filter_values))
      for(fi in 1:length(filter_values)) {
        if(fi == 1 & filter_values[1] == minimal_value){
          if(filter_NA) {	name_list[1] <- "non-missing"
          } else {    		name_list[1] <- "all" }
        } else {
          if(filter_values[fi] < 1) {	name_list[fi] <- paste0(var_name, sign, " ", filter_values[fi])
          } else	{	      			    	name_list[fi] <- paste0(var_name, sign, " ", filter_values[fi], "/N") }
        }
      }
      return(name_list)
    }
    
    QQfilter_data <- function(input_var, input_N, filter_values, filter_NA = TRUE, minimal_value = 0, two_sided = FALSE, subset = !logical(length = length(input_var))) {
      filter_list <- matrix(data = FALSE, nrow = length(input_var), ncol = length(filter_values))
      if(filter_NA) {
        var_NA <- is.na(input_var)
        if(two_sided) {
          for(fi in 1:length(filter_values)) {
            if(filter_values[fi] == minimal_value) { filter_list[ , fi] <- var_NA & subset
            } else if(filter_values[fi] < 1) {       filter_list[ , fi] <- (var_NA | input_var < filter_values[fi]	| input_var > 1 - filter_values[fi]) & subset
            } else {				                    		 filter_list[ , fi] <- (var_NA | input_var < filter_values[fi]/input_N | input_var > 1 - filter_values[fi]/input_N | is.na(input_N) ) & subset }
          }
        } else { #not two_sided
          for(fi in 1:length(filter_values)) {
            if(filter_values[fi] == minimal_value) { filter_list[ , fi] <- var_NA & subset
            } else if(filter_values[fi] < 1) { 	     filter_list[ , fi] <- (var_NA | input_var < filter_values[fi]) & subset
            } else {						                  	 filter_list[ , fi] <- (var_NA | input_var < filter_values[fi]/input_N | is.na(input_N) ) & subset }
          }
        } # end two-sided if
      } else {	# when filter_NA is FALSE
        var_nonNA <- !is.na(input_var)
        if(two_sided) {
          for(fi in 1:length(filter_values)) {
            if(filter_values[fi] != minimal_value) { # if it is the minimal value, the filter is ignored
              if(filter_values[fi] < 1) {            filter_list[ , fi] <- var_nonNA & (input_var < filter_values[fi] | input_var > 1 - filter_values[fi]) & subset
              } else {					                     filter_list[ , fi] <- var_nonNA & (input_var < filter_values[fi]/input_N | input_var > 1 - filter_values[fi]/input_N ) & subset & !is.na(input_N) }
            }
          }
        } else { #not two_sided
          for(fi in 1:length(filter_values)) {
            if(filter_values[fi] != minimal_value) { # if it is the minimal value, the filter is ignored
              if(filter_values[fi] < 1) {            filter_list[ , fi] <- var_nonNA & (input_var < filter_values[fi] ) & subset
              } else {					                     filter_list[ , fi] <- var_nonNA & (input_var < filter_values[fi]/input_N ) & subset & !is.na(input_N) }
            }
          }
        } # end two-sided if
      }
      return(filter_list)
    }
    
    if(ignore_impstatus & !is.null(c(HWEfilter_values, calfilter_values, impfilter_values))) {
      no_subset <- !logical(length = dataN) }
  }
  
  # Filter for allele frequency
  if(is.null(FRQfilter_values)){
    useFRQfilter <- FALSE
    FRQfilter_out <- "no MAF filter applied"
    FRQfilter_N	 <- NULL
  } else {
    FRQfilter_values	<- QQfilter_sort(FRQfilter_values)
    FRQfilter		<- QQfilter_data(dataset$EFF_ALL_FREQ, input_N = dataset$N_TOTAL, filter_values = FRQfilter_values, filter_NA = FRQfilter_NA, two_sided = TRUE)
    FRQfilter_N		<- colSums(FRQfilter)
    FRQfilter_count	<- length(FRQfilter_values)
    useFRQfilter	<- FRQfilter_N > c(0, FRQfilter_N[1:(FRQfilter_count - 1)])
    
    if(any(FRQfilter_values >= 1) & which.max(FRQfilter_values) < FRQfilter_count) {
      if(FRQfilter_NA) { useFRQfilter[which.max(FRQfilter_values) + 1] <- FRQfilter_N[which.max(FRQfilter_values) + 1] > FRQfilter_N[1]
      } else {		 useFRQfilter[which.max(FRQfilter_values) + 1] <- FRQfilter_N[which.max(FRQfilter_values) + 1] > 0 }
    }
    
    FRQfilter_names	<- paste0(QQfilter_name("", filter_values = FRQfilter_values, filter_NA = FRQfilter_NA),
                              " (", round( ( 1 - FRQfilter_N/dataN) * 100), "%)")
    FRQfilter_out	<- QQfilter_name("MAF ", filter_values = FRQfilter_values, filter_NA = FRQfilter_NA, sign = "<")
  }
  
  # Filter for HWE p-values
  if(is.null(HWEfilter_values)){
    useHWEfilter <- FALSE
    HWEfilter_out <- "no HWE p-value filter applied"
    HWEfilter_N	 <- NULL
  } else {
    HWEfilter_values	<- QQfilter_sort(HWEfilter_values)
    HWEfilter		<- QQfilter_data(dataset$HWE_PVAL, input_N = dataset$N_TOTAL, filter_values = HWEfilter_values, filter_NA = HWEfilter_NA, subset = if(ignore_impstatus) no_subset else genoset)
    HWEfilter_N		<- colSums(HWEfilter)
    HWEfilter_count	<- length(HWEfilter_values)
    useHWEfilter	<- HWEfilter_N > c(0, HWEfilter_N[1:(HWEfilter_count - 1)])
    
    if(any(HWEfilter_values >= 1) & which.max(HWEfilter_values) < HWEfilter_count) {
      if(HWEfilter_NA) { useHWEfilter[which.max(HWEfilter_values) + 1] <- HWEfilter_N[which.max(HWEfilter_values) + 1] > HWEfilter_N[1]
      } else {		 useHWEfilter[which.max(HWEfilter_values) + 1] <- HWEfilter_N[which.max(HWEfilter_values) + 1] > 0 }
    }
    
    HWEfilter_names	<- paste0(QQfilter_name("", filter_values = HWEfilter_values, filter_NA = HWEfilter_NA),
                              " (", round( ( 1 - HWEfilter_N/dataN) * 100), "%)")
    HWEfilter_out	<- QQfilter_name("HWE p-value ", filter_values = HWEfilter_values, filter_NA = HWEfilter_NA, sign = "<")
  }
  
  # Filter for call rate
  if(is.null(calfilter_values)){
    useCalfilter <- FALSE
    calfilter_out <- "no call rate filter applied"
    calfilter_N	 <- NULL
  } else {
    calfilter_values	<- QQfilter_sort(calfilter_values)
    calfilter		<- QQfilter_data(dataset$CALLRATE, input_N = dataset$N_TOTAL, filter_values = calfilter_values, filter_NA = calfilter_NA, subset = if(ignore_impstatus) no_subset else genoset)
    calfilter_N		<- colSums(calfilter)
    calfilter_count	<- length(calfilter_values)
    useCalfilter	<- calfilter_N > c(0, calfilter_N[1:(calfilter_count - 1)])
    
    if(any(calfilter_values >= 1) & which.max(calfilter_values) < calfilter_count) {
      if(calfilter_NA) { useCalfilter[which.max(calfilter_values) + 1] <- calfilter_N[which.max(calfilter_values) + 1] > calfilter_N[1]
      } else {		 useCalfilter[which.max(calfilter_values) + 1] <- calfilter_N[which.max(calfilter_values) + 1] > 0 }
    }
    
    calfilter_names	<- paste0(QQfilter_name("", filter_values = calfilter_values, filter_NA = calfilter_NA),
                              " (", round( ( 1 - calfilter_N/dataN) * 100), "%)")
    calfilter_out	<- QQfilter_name("call rate ", filter_values = calfilter_values, filter_NA = calfilter_NA, sign = "<")
  }
  
  # Filter for imputation quality
  if(is.null(impfilter_values)){
    useImpfilter <- FALSE
    impfilter_out <- "no imputation-quality filter applied"
    impfilter_N	 <- NULL
  } else {
    impfilter_values	<- QQfilter_sort(impfilter_values, minimal_value = impfilter_min)
    impfilter		<- QQfilter_data(dataset$IMP_QUALITY, input_N = dataset$N_TOTAL, filter_values = impfilter_values, filter_NA = impfilter_NA, minimal_value = impfilter_min, subset = if(ignore_impstatus) no_subset else imp_set)
    impfilter_N		<- colSums(impfilter)
    impfilter_count	<- length(impfilter_values)
    useImpfilter	<- impfilter_N > c(0, impfilter_N[1:(impfilter_count - 1)])
    
    if(any(impfilter_values >= 1) & which.max(impfilter_values) < impfilter_count) {
      if(impfilter_NA) { useImpfilter[which.max(impfilter_values) + 1] <- impfilter_N[which.max(impfilter_values) + 1] > impfilter_N[1]
      } else {		 useImpfilter[which.max(impfilter_values) + 1] <- impfilter_N[which.max(impfilter_values) + 1] > 0 }
    }
    
    impfilter_names	<- paste0(QQfilter_name("", filter_values = impfilter_values, filter_NA = impfilter_NA, minimal_value = impfilter_min),
                              " (", round( ( 1 - impfilter_N/dataN) * 100), "%)")
    impfilter_out	<- QQfilter_name("Imputation Q ", filter_values = impfilter_values, filter_NA = impfilter_NA, minimal_value = impfilter_min, sign = "<")
  }
  
  
  # Calculating lambda
  # First the script creates a short-list by throwing out NA's (!) to calculate the lambda's,
  #	and a p-cutoff filter is created (!). The script checks whether there are still more than
  #	10 datapoints available at both exclamation marks. If yes, the script uses the shortlist
  #	to generate the unfiltered expected and observed QQ plots; and the long list plus the 
  #	addQQplot function to generate the filtered plots.
  #	WARNING: note that the content of QQ_obs_short changes: it starts out as the p-values
  #	minus NA (*), then -log10 values (#), then -log10 values minus p > 0.05 (***). Lambda needs
  #	to be calculated at the (*). QQ_exp is calculated at (#), but (*) works too; the contents
  #	of the vector used to calculate QQ_exp do not matter, only its length. QQ_exp changes
  #	as well: a copy without > 0.05 is made at (***: QQ_exp_short). The copy is used for
  #	the graphs, but QQ_exp is necessary for the probability clouds. If so, it's shortened
  #	to 1000 points to prevent the polygon functions from being overloaded.
  
  QQ_obs_p	 <- dataset$PVALUE
  QQ_obs_short <- subset(QQ_obs_p, !is.na(QQ_obs_p))
  QQ_obs_N	 <- length(QQ_obs_short)	# N of non-missing p's: the N of all p's is dataN
  
  if(QQ_obs_N > 10L) {
    lambda <- median(qchisq(QQ_obs_short, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
    
    if(no_impstatus) {
      lambda_geno <- NA
      lambda_imp  <- NA
    } else {
      if(sum(genoset & !is.na(QQ_obs_p)) > 10L) { lambda_geno <- median(qchisq(subset(QQ_obs_p, genoset & !is.na(QQ_obs_p)), df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
      } else { lambda_geno <- NA }
      if(sum(imp_set & !is.na(QQ_obs_p)) > 10L) { lambda_imp	<- median(qchisq(subset(QQ_obs_p,	imp_set & !is.na(QQ_obs_p)), df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
      } else { lambda_imp	<- NA }
    }
    
    QQ_obs_short<- sort(-log10(QQ_obs_short))
    QQ_incl	<- QQ_obs_short >= -log10(plot_cutoff_p)
  } else {
    lambda	<- NA
    lambda_geno <- NA
    lambda_imp	<- NA
    QQ_incl	<- 0
  }
  
  # Stage 2: creating QQ plots
  if(plot_QQ & sum(QQ_incl) > 10L ) {
    print(" - creating QQ plots", quote = FALSE)
    flush.console()
    
    # creating functions for QQ plots
    calculateQQ <- function(p_length) { return(sort(-log10(ppoints(p_length)))) }
    # The ppoints function (and hence the calculate QQ function)
    # accepts both the array of p-values, or just its length
    
    addQQplot <- function(p_values, cutoff_p, input_filter, input_color) {
      p_values <- p_values[!input_filter & !is.na(p_values)] 
      n <- length(p_values)
      if(n > 0L) {
        p_values <- sort(-log10(p_values))
        incl <- p_values >= -log10(cutoff_p)
        if(incl[n]) points(calculateQQ(n)[incl], p_values[incl], pch = 1, col = input_color)
      }
      return(invisible())
    }
    
    QQ_exp	<- calculateQQ(QQ_obs_N)
    QQ_exp_short<- QQ_exp[QQ_incl]
    QQ_obs_short<- QQ_obs_short[QQ_incl]
    QQ_exp_min	<- QQ_exp_short[1]
    QQ_exp_max	<- QQ_exp_short[length(QQ_exp_short)]
    QQ_obs_min	<- QQ_obs_short[1]
    QQ_obs_max	<- QQ_obs_short[length(QQ_obs_short)]
    
    if(plot_QQ_bands) {
      temp <- (1:QQ_obs_N)
      i1000 <- c(1, (1:1000) * floor(QQ_obs_N / 1000), QQ_obs_N)
      QQ_band_upper <- sort(-log10(qbeta( 1 - 0.05 / 2, temp, QQ_obs_N - temp + 1 ) ) )[i1000]
      QQ_band_lower <- sort(-log10(qbeta(		 0.05 / 2, temp, QQ_obs_N - temp + 1 ) ) )[i1000]
      QQ_exp <- QQ_exp[i1000]
    }
    
    # If no filter values have been specified, the FRQ plot is printed by default (with a different header)
    QQplot_N <- sum(!is.null(FRQfilter_values), !is.null(HWEfilter_values), !is.null(calfilter_values), !is.null(impfilter_values))
    if(QQplot_N == 0L) QQplot_N <- 1L

    png(paste0(save_dir, "/", save_name, "_graph_QQ.png"),
        width = if(QQplot_N == 4L) 1440 else QQplot_N * 720,
        height= if(QQplot_N == 4L) 1440 else 720, res = 144)
    par(mfrow = if(QQplot_N == 4L) c(2,2) else c(1, QQplot_N ))
    
    
    # allele frequency QQ plot
    if(!is.null(FRQfilter_values) | is.null(c(HWEfilter_values, calfilter_values, impfilter_values)) ) {
      plot(c(QQ_exp_min, QQ_exp_max), c(QQ_obs_min, QQ_obs_max), xlim = c(0, QQ_exp_max), ylim = c(0, QQ_obs_max),
           main = if(is.null(FRQfilter_values)) { "QQ plot" } else { "QQ plot - allele frequency" },
           xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)",
           pch = 1, col = "black")
      mtext("Imputed & genotyped SNPs", cex = 0.8, line = 0.5, col="blue")
      if(plot_QQ_bands) {
        polygon( c(QQ_exp, rev(QQ_exp)), c( QQ_band_upper, rev(QQ_exp)), col="grey", border = NA )
        polygon( c(QQ_exp, rev(QQ_exp)), c( QQ_band_lower, rev(QQ_exp)), col="grey", border = NA )
      }
      points(QQ_exp_short, QQ_obs_short, pch = 1, col = "black")
      if(any(useFRQfilter)) {
        for(ip in 1:FRQfilter_count) {
          if(useFRQfilter[ip]) addQQplot(QQ_obs_p, plot_cutoff_p, FRQfilter[ ,ip], QQ_colors[ip])
        }
      }
      if(!is.null(FRQfilter_values)) {
        legend(0, 0.94*QQ_obs_max, c("All", c(FRQfilter_names)[useFRQfilter] ),
               pch = 1, col = c("black", QQ_colors[useFRQfilter]))
      }
      abline(0,1)
      text(0, 0.98 * QQ_obs_max, substitute(paste(lambda[GC], "=", x), list(x=round(lambda, digits = 3))), pos = 4, col = ifelse(lambda > 1.1, "red", "black"))
    }
    
    # HWE p-value QQ plot
    if(!is.null(HWEfilter_values)) {
      plot(c(QQ_exp_min, QQ_exp_max), c(QQ_obs_min, QQ_obs_max), xlim = c(0, QQ_exp_max), ylim = c(0, QQ_obs_max),
           main = "QQ plot - HWE p-value", xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)",
           pch = 1, col = "black")
      if(ignore_impstatus) { mtext("Imputed & genotyped SNPs", cex = 0.8, line = 0.5, col="red")
      } else {         mtext("Genotyped SNPs only", cex = 0.8, line = 0.5, col="blue") }
      if(plot_QQ_bands) {
        polygon( c(QQ_exp, rev(QQ_exp)), c( QQ_band_upper, rev(QQ_exp)), col="grey", border = NA )
        polygon( c(QQ_exp, rev(QQ_exp)), c( QQ_band_lower, rev(QQ_exp)), col="grey", border = NA )
      }
      points(QQ_exp_short, QQ_obs_short, pch = 1, col = "black")
      if(any(useHWEfilter)) {
        for(ip in 1:HWEfilter_count) {
          if(useHWEfilter[ip]) addQQplot(QQ_obs_p, plot_cutoff_p, HWEfilter[ ,ip], QQ_colors[ip])
        }
      }
      legend(0, 0.94*QQ_obs_max, c("All", c(HWEfilter_names)[useHWEfilter] ),
             pch = 1, col = c("black", QQ_colors[useHWEfilter]))
      abline(0,1)
      if(is.null(FRQfilter_values)) {
        text(0, 0.98 * QQ_obs_max, substitute(paste(lambda[GC], "=", x), list(x=round(lambda, digits = 3))), pos = 4, col = ifelse(lambda > 1.1, "red", "black"))
      }
    }
    
    
    # call rate QQ plot
    if(!is.null(calfilter_values)) {
      plot(c(QQ_exp_min, QQ_exp_max), c(QQ_obs_min, QQ_obs_max), xlim = c(0, QQ_exp_max), ylim = c(0, QQ_obs_max),
           main = "QQ plot - call rates", xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)",
           pch = 1, col = "black")
      if(ignore_impstatus) { mtext("Imputed & genotyped SNPs", cex = 0.8, line = 0.5, col="red")
      } else {         mtext("Genotyped SNPs only", cex = 0.8, line = 0.5, col="blue") }
      if(plot_QQ_bands) {
        polygon( c(QQ_exp, rev(QQ_exp)), c( QQ_band_upper, rev(QQ_exp)), col="grey", border = NA )
        polygon( c(QQ_exp, rev(QQ_exp)), c( QQ_band_lower, rev(QQ_exp)), col="grey", border = NA )
      }
      points(QQ_exp_short, QQ_obs_short, pch = 1, col = "black")
      if(any(useCalfilter)) {
        for(ip in 1:calfilter_count) {
          if(useCalfilter[ip]) addQQplot(QQ_obs_p, plot_cutoff_p, calfilter[ ,ip], QQ_colors[ip])
        }
      }
      legend(0, 0.94*QQ_obs_max, c("All", c(calfilter_names)[useCalfilter] ),
             pch = 1, col = c("black", QQ_colors[useCalfilter]))
      abline(0,1)
      if(is.null(FRQfilter_values) & is.null(HWEfilter_values)) {
        text(0, 0.98 * QQ_obs_max, substitute(paste(lambda[GC], "=", x), list(x=round(lambda, digits = 3))), pos = 4, col = ifelse(lambda > 1.1, "red", "black"))
      }
    }

    
    # imputation quality QQ plot
    if(!is.null(impfilter_values)) {
      plot(c(QQ_exp_min, QQ_exp_max), c(QQ_obs_min, QQ_obs_max), xlim = c(0, QQ_exp_max), ylim = c(0, QQ_obs_max),
           main = "QQ plot - imputation quality", xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)",
           pch = 1, col = "black")
      if(ignore_impstatus) { mtext("Imputed & genotyped SNPs", cex = 0.8, line = 0.5, col="red")
      } else {         mtext("Imputed SNPs only", cex = 0.8, line = 0.5, col="blue") }
      if(plot_QQ_bands) {
        polygon( c(QQ_exp, rev(QQ_exp)), c( QQ_band_upper, rev(QQ_exp)), col="grey", border = NA )
        polygon( c(QQ_exp, rev(QQ_exp)), c( QQ_band_lower, rev(QQ_exp)), col="grey", border = NA )
      }
      points(QQ_exp_short, QQ_obs_short, pch = 1, col = "black")
      if(any(useImpfilter)) {
        for(ip in 1:impfilter_count) {
          if(useImpfilter[ip]) addQQplot(QQ_obs_p, plot_cutoff_p, impfilter[ , ip], QQ_colors[ip])
        }
      }
      legend(0, 0.94*QQ_obs_max, c("All", c(impfilter_names)[useImpfilter] ),
             pch = 1, col = c("black", QQ_colors[useImpfilter]))
      abline(0,1)
      if(QQplot_N == 1) {
        text(0, 0.98 * QQ_obs_max, substitute(paste(lambda[GC], "=", x), list(x=round(lambda, digits = 3))), pos = 4, col = ifelse(lambda > 1.1, "red", "black"))
      } 
    }
    title(sub = save_name, cex.sub = 1.3)
    dev.off()
  } else {
    print(" - - creating QQ plots (SKIPPED)", quote = FALSE)
    if(plot_QQ) {
      if(use_log) save_log(phaseL = 4L, checkL = "QQ plots", typeL = "insufficient data", SNPL = sum(QQ_incl), allSNPs = dataN, actionL = "Test skipped", noteL = "Insufficient significant p-values to create QQ plots", fileL = paste(save_dir, save_name, sep = "/"))
      print(paste0(" - - insufficient significant p-values (", sum(QQ_incl), " SNP with p =< ", plot_cutoff_p, ") to create QQ plots"), quote = FALSE)
    }
  }
  
  # Stage 3: creating Manhattan plot
  if(plot_Man & sum(QQ_incl) > 10L ) {
    print(" - creating Manhattan plot", quote = FALSE)
    flush.console()
    
    manfilter <- HQ_filter(data = dataset, ignore_impstatus = ignore_impstatus,
                           FRQ_val = manfilter_FRQ, HWE_val = manfilter_HWE, cal_val = manfilter_cal, imp_val = manfilter_imp,
                           FRQ_NA = FRQfilter_NA, HWE_NA = HWEfilter_NA, cal_NA = calfilter_NA, imp_NA = impfilter_NA)
    manfilter_N <- dataN - sum(manfilter)
    if(manfilter_N == dataN) {
      manhattanN <- 0L
    } else {
      man_set <- dataset[manfilter & !is.na(dataset$POSITION) & !is.na(dataset$CHR) & !is.na(dataset$PVALUE) & dataset$PVALUE <= plot_cutoff_p, c("MARKER", "CHR", "POSITION", "PVALUE")[c(plot_names, TRUE, TRUE, TRUE)]]
      
      if(!is.numeric(man_set$CHR)) {
        man_set$CHR <- toupper(man_set$CHR)
        man_set$CHR[man_set$CHR == "X"] <- 23L
        man_set$CHR[man_set$CHR == "Y"] <- 24L
        man_set$CHR[man_set$CHR == "XY"] <- 25L
        man_set$CHR[man_set$CHR %in% c("M", "MT")] <- 26L
        
        man_set$CHR <- as.integer(man_set$CHR)
        if(any(is.na(man_set$CHR))) {
          print(" - - Cannot translate chromosome values - entries excluded from Manhattan plot", quote = FALSE)
          man_set <- man_set[!is.na(man_set$CHR), ]
        }
      }
      
      if(any(!man_set$CHR %in% 1:26)) {
        print(" - - chromosome values outside normal range - entries excluded from Manhattan plot", quote = FALSE)
        man_set <- man_set[man_set$CHR %in% 1:26, ]				
      }
      manhattanN <- nrow(man_set)
    }
    
    if(manhattanN > 9L) { 		# testing if there are sufficient p-values < 0.05
      chr_size <- data.frame(chromosome = 1:27, #	chr2				chr3			chr4				chr5				chr6				chr7				chr8				chr9			chr10				chr11			 chr12			chr13				chr14				chr15				chr16			chr17			chr18				chr19			chr20				chr21			chr22			X23					Y24		XY25,M26, 27 = end of M
                             size = c(247249719, 242951149, 199501827, 191273063, 180857866,	170899992,	158821424,	146274826,	140273252,	135374737,	134452384,	132349534,	114142980,	106368585,	100338915,	 88827254,	 78774742,	 76117153,	 63811651,	 62435964,	 46944323,	 49691432,	154913754,	 57772954,	0,	0,	0 ),
                             start= c(	 500000, 248249719, 491700868, 691702695, 883475758, 1064833624, 1236233616, 1395555040, 1542329866, 1683103118, 1818977855, 1953930239, 2086779773, 2201422753, 2308291338, 2409130253, 2498457507, 2577732249, 2654349402, 2718661053, 2781597017, 2829041340, 2879232772, 3034646526, NA, NA, NA))
      
      new_pos <- integer(length = manhattanN)
      for (mi in 1:24 ) new_pos <- ifelse(man_set$CHR==chr_size$chromosome[mi], chr_size[mi, 3] + man_set$POSITION, new_pos)
      
      use_Y	<- any(man_set$CHR == 24L)
      use_XY <- any(man_set$CHR == 25L)
      use_M	<- any(man_set$CHR == 26L)
      man_label <- c(1:22, "X")
      at_label	<- chr_size$start[2:24] - 250000 - ( chr_size$start[2:24] - chr_size$start[1:23] ) / 2
      #		at_label <-	 chr_size$start[1:23] + (chr_size$start[2:24] - 500000 - chr_size$start[1:23])/2
      
      if(use_Y) {
        man_label <- c(man_label, "Y")
        at_label <- c(at_label, chr_size$start[25] - 250000 - ( chr_size$start[25] - chr_size$start[24] ) / 2)
        chr_size$start[25] <- chr_size$start[24] + chr_size$size[24] + 500000
      } else { chr_size$start[25] <- chr_size$start[24] }
      
      if(use_XY) {
        man_label <- c(man_label, "XY")
        at_label <- c(at_label, chr_size$start[26] - 250000 - ( chr_size$start[26] - chr_size$start[25] ) / 2)
        chr_size$size[25] <- max(man_set$POSITION[man_set$CHR == 25])
        chr_size$start[26] <- chr_size$start[25] + chr_size$size[25] + 500000
      } else { chr_size$start[26] <- chr_size$start[25] }
      
      if(use_M ) {
        man_label <- c(man_label, "M")
        at_label <- c(at_label, chr_size$start[27] - 250000 - ( chr_size$start[27] - chr_size$start[26] ) / 2)
        chr_size$size[26] <- max(man_set$POSITION[man_set$CHR == 26])
        chr_size$start[27] <- chr_size$start[26] + chr_size$size[26] + 500000
      } else { chr_size$start[27] <- chr_size$start[26] }
      
      manMax <- ceiling(-log10(min(man_set$PVALUE)))
      if(manMax < 10L) manMax <- 10L
      
      png(paste0(save_dir, "/", save_name, "_graph_M.png"),
           width = 960, height = 480)
      par(mfrow = c(1,1), mgp = c(2.5, 0.9, 0))
      palette(c("red", "green3", "blue", "cyan"))
      plot(new_pos, -log10(man_set$PVALUE), pch = 20, xaxs = "i", xaxt = "n", ylim = c(1, manMax),
           xlim = c(0, chr_size$start[27]), col = man_set$CHR, cex.lab = 1.8, cex.axis = 1.4,
           xlab = "Chromosome", ylab = "Observed -log10(p-value)", main = "Manhattan plot", sub = save_name, cex.sub = 1.3, col.sub = "red")
      abline(v = chr_size[2:26, 3] - 250000, lty = 2)
      abline(h = -log10(5e-8), lty = 3, col="red")
      axis(1, at = at_label, labels = man_label, cex.axis = 1.5)
      dev.off()
    } else {
      if(use_log) save_log(phaseL = 4L, checkL = "Manhattan plot", typeL = "insufficient data", SNPL = manhattanN, allSNPs = dataN, actionL = "-", noteL = "Insufficient signif., unfiltered SNPs with known locations to create Manhattan plot", fileL = paste(save_dir, save_name, sep = "/"))
      print(paste0(" - - insufficient significant, unfiltered SNPs (N = ", manhattanN, ") with known location to create Manhattan plot"), quote = FALSE)
    }
  } else {
    print(" - - creating Manhattan plot (SKIPPED)", quote = FALSE)
    if(plot_Man) {
      if(use_log) { save_log(phaseL = 4L, checkL = "Manhattan plot", typeL = "insuf. data", SNPL = sum(QQ_incl), allSNPs = dataN, actionL = "Test skipped", noteL = "Insufficient significant p-values to create Manhattan plot", fileL = paste(save_dir, save_name, sep = "/"))
      } else { print(paste0(" - - insufficient significant p-values (", sum(QQ_incl), " SNP with p < ", plot_cutoff_p, ") to create Manhattan plot"), quote = FALSE) }
    }
    manfilter_N <- NA
    manhattanN	<- NA
  }
  flush.console()
  return(list(lambda = c(lambda, lambda_geno, lambda_imp), ignore_impstatus = ignore_impstatus,
              FRQfilter_names = FRQfilter_out, HWEfilter_names = HWEfilter_out, calfilter_names = calfilter_out, impfilter_names = impfilter_out,
              FRQfilter_N = FRQfilter_N, HWEfilter_N = HWEfilter_N, calfilter_N = calfilter_N, impfilter_N = impfilter_N,
              Manfilter_N = manfilter_N))
}
