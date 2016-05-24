
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("lodGWAS library, version 1.0-7")
  packageStartupMessage("")
  packageStartupMessage("The user manual can be found in")
  packageStartupMessage(system.file("doc", package = "lodGWAS"))
}

# FUNCTION INFORMATION COMPUTES THE IMPUTATION QUALITY SCORE OF THE SNP
imputation_quality <- function(dosage) {
  n2 <- (sum(dosage^2) - sum(dosage)) / 2
  n1 <- sum(dosage) - 2*n2
  n0 <- length(dosage) - n1 -n2
  
  f0 <- (n0+n1/2) / (n0+n1+n2)
  return(if(f0>0 & f0<1) { sd(dosage)^2 / (2*f0*(1-f0)) } else 0)
} # END FUNCTION INFORMATION


# FUNCTION SURVIVAL PERFORMS A LEFT AND RIGHT CENSORED SURVIVAL ANALYSIS 
# TO DEAL WITH VALUES OUTSIDE THE DETECTION INTERVAL
survivalf <- function(trait, outsideLOD, basic_model=NULL, dist="gaussian",
                      SNP, data, filterpar=NULL, filterval=NULL) {
  
  # Checking & defining data
  checkTrait <- function(values, name){
    if(!is.vector(values)) stop(paste0("'", name, "' is not a vector"))
    if(length(values) == 0L) stop(paste0("The variable '", name ,"' has length 0"))
    if(!is.numeric(values)) stop(paste0("'", name, "' is not a numeric variable"))
    return(invisible(NULL))
  }
  
  checkOutsideLOD <- function(values, traitvalues){
    if(!is.vector(values)) stop("'outsideLOD' must be a vector")
    if(length(values) == 0L) stop("The variable 'outsideLOD' has length 0")
    if(!is.integer(values)) stop("'outsideLOD' must be an integer variable")
    if(length(values) != length(traitvalues)) stop("Variables 'trait' and 'outsideLOD' are of unequal length")
    if(any(is.na(values) & !is.na(traitvalues))) stop("Missing values in 'outsideLOD'")
    if(any(values < 0L | values > 2L, na.rm = TRUE)) stop("'outsideLOD' contains invalid (i.e. not 0, 1 or 2) values")
    return(invisible(NULL))
  }
  
  if(missing(data)) {
    if(missing(trait)) stop("Trait variable not found")
    if(missing(outsideLOD)) stop("outsideLOD variable not found")
    if(missing(SNP)) stop("SNP variable not found")
    
    checkTrait(trait, "trait")
    checkOutsideLOD(outsideLOD, trait)
    checkTrait(SNP, "SNP")
    if(length(trait) != length(SNP)) stop("Variables 'trait' and 'SNP' are of unequal length")
  
    data <- data.frame(trait = trait, outsideLOD = outsideLOD, SNP = SNP)
    rm(trait, outsideLOD, SNP)
  } else {
    if(missing(trait)) {
      if(!"trait" %in% colnames(data)) stop("No variable 'trait' defined or present in data")
      checkTrait(data$trait, "trait")
    } else {
      checkTrait(trait, "trait")
      if(length(trait) != nrow(data)) stop("Variables 'trait' and 'data' are of unequal length")
      data$trait <- trait
      rm(trait)
    }
    
    if(missing(outsideLOD)) {
      if(!"outsideLOD" %in% colnames(data)) stop("No variable 'outsideLOD' defined or present in data")
      checkOutsideLOD(data$outsideLOD, data$trait)
    } else {
      checkOutsideLOD(outsideLOD, data$trait)
      data$outsideLOD <- outsideLOD
      rm(outsideLOD)
    }
    
    if(missing(SNP)) {
      if(!"SNP" %in% colnames(data)) stop("No variable 'SNP' defined or present in data")
      checkTrait(data$SNP, "SNP")
    } else {
      checkTrait(SNP, "SNP")
      if(length(SNP) != nrow(data)) stop("Variables 'SNP' and 'data' are of unequal length")
      data$SNP <- SNP
      rm(SNP)
    }
  }
  
  # Running filter
  if (!is.null(filterpar)) {
    if (is.null(filterval)) {
      cat("Warning: you specified a filter variable, but no filter value.\n         None of the observations will be removed.\n")
    } else {
      if (is.na(filterval)) { 
        cat("Warning: you specified a filter variable, but no filter value.\n         None of the observations will be removed.\n")
      } else {
        NNApre <- sum(is.na(data$trait))
        if(length(filterpar) == 1){
          if(!filterpar %in% colnames(data)) stop("Cannot find variable 'filterpar' in 'data'")
          data$trait[which(data[ , filterpar] != filterval)] <- NA
          
        } else {
          if(length(filterpar) != nrow(data)) stop("Variable 'filterpar' is not of equal length to data")
          data$trait[which(filterpar != filterval)] <- NA
        }
        if(all(is.na(data$trait))) stop("After filtering, no values are left.")
        NNApost <- sum(is.na(data$trait))
        print(paste("Filtering removed", NNApost - NNApre, "values from the dataset"), quote = FALSE)        
      } } }
  
  # start of actual analysis
  # result consists of 7 items: 1 = N-total; 2 = N-valid; 3 = coded allele frequency;
  #                             4 = beta; 5 = SE; 6 = p; 7 = imp info
  result <- rep(NA,7)
  
  AF <- sum(data$SNP)/(2*length(data$SNP))
  info <- imputation_quality(data$SNP)
  result[3] <- AF
  result[7] <- info
  if (AF>0.001 & AF<0.999 & info>0.01) {
    S <- Surv(data$trait, data$trait, data$outsideLOD, type="interval")
    model <- if(is.null(basic_model)) "S~SNP" else paste0("S~", basic_model, "+SNP")
    options(warn=-1) # temporarily disables warnings
    m <- psm(as.formula(model), data=data, dist=dist)
    options(warn=0) # reenables warnings
    result[4] <- m$coefficients["SNP"]
    result[5] <- sqrt(m$var["SNP","SNP"])
    if (m$coefficients["SNP"]!=0   &   sqrt(m$var["SNP", "SNP"])!=0) {
      result[6] <- 2 * pnorm(-abs(m$coefficients["SNP"]) / sqrt(m$var["SNP","SNP"]))
    }
    result[1] <- m$stats[1]
    
    filterLOD <- which(data$outsideLOD==1L & !is.na(data$trait))
    data <- data[filterLOD,]
    
    # the commands below calculate N. S needs to be recalculated because it is part
    # of "model". Removing the command causes error messages
    S <- Surv(data$trait, data$outsideLOD)
    m <- psm(as.formula(model), data=data, dist="gaussian")
    result[2] <- m$stats[1]
  }
  
  return(result)
} # END FUNCTION SURVIVAL


lod_QC <- function(phenofile, pheno_name,
                  filedirectory = getwd(),
                  outputfile = "lodQC",
                  lower_limit = NA, upper_limit = NA,
                  stop_if_error = FALSE, reprint_warnings = FALSE){
  #if(is.na(lower_limit) & is.na(upper_limit)) stop("No limits specified")
  
  stopifnot(is.character(filedirectory), is.character(filedirectory), length(filedirectory) == 1L)
  if(!file.exists(filedirectory)) stop("Cannot find 'filedirectory'")
  
  if(missing(pheno_name)) pheno_name <- "trait"

  stopifnot(is.vector(outputfile), is.character(outputfile), length(outputfile) == 1L, nchar(outputfile) > 0L)
  
  ### loading phenofile
  if(!is.data.frame(phenofile) & !is.matrix(phenofile)) {
    if(!is.character(phenofile) | length(phenofile) != 1L) stop("'phenofile' is neither a dataset nor a filename")
    if(!file.exists(paste0(filedirectory, "/", phenofile))) stop(paste0("Cannot find file '", phenofile, "' in the working directory"))
    phenofile <- read.table(file = paste0(filedirectory, "/", phenofile), stringsAsFactors = FALSE, header = TRUE, comment.char = "")
  }
  
  if(!all(c(pheno_name, "outsideLOD") %in% colnames(phenofile))) stop(paste0("Cannot find columns 'outsideLOD' and/or '", pheno_name, "' in 'phenofile'"))
  if(!is.numeric(phenofile[ , pheno_name])) stop(paste0("Column '", pheno_name, "' does not contain numeric values"))
  if(!is.integer(phenofile$outsideLOD)) stop("Column 'outsideLOD' does not contain integer values")
  
  ### collecting sample sizes BEFORE removal of missing phenos (N0)
  N0_all <- nrow(phenofile)
  N0_NA <- sum(is.na(phenofile$outsideLOD))
  N0_lod0 <- sum(phenofile$outsideLOD == 0, na.rm = TRUE)
  N0_lod1 <- sum(phenofile$outsideLOD == 1, na.rm = TRUE)
  N0_lod2 <- sum(phenofile$outsideLOD == 2, na.rm = TRUE)
  
  ### removing missing phenos (N1)
  phenofile <- phenofile[!is.na(phenofile[ , pheno_name]), c(pheno_name, "outsideLOD") ]
  colnames(phenofile)[1] <- "trait"
  N1_all <- nrow(phenofile)
  if(N1_all == 0L) stop("No non-missing phenotypes present")
  if(any(is.na(phenofile$outsideLOD))) stop("Column 'outsideLOD' contains missing values")
  if(any(phenofile$outsideLOD > 2L | phenofile$outsideLOD < 0L)) stop("Column 'outsideLOD' can only contain values 0, 1 or 2")
  N1_lod1 <- sum(phenofile$outsideLOD == 1L)
  if(N1_lod1 == 0L) stop("Column 'outsideLOD' marks all values as outside of LOD")
  
  
  ### checking lower limit
  
  N1_lod2 <- sum(phenofile$outsideLOD == 2L)
  if(is.na(lower_limit)){
    N1_llim_lod1 <- 0L
    N1_llim_lod2 <- 0L
    N1_llim4 <- 0L
  } else {
    if(!is.na(upper_limit) & lower_limit >= upper_limit) stop("'upper_limit' does not exceed 'lower_limit'")
    # llim1: are there pheno entries < LL
    # llim2: are there phenos < LL not marked as LOD2
    # llim3: are there lod2 entries that exceed LL
    # llim4: are there phenoentries that equal LL, yet are LOD = 0
    # llim_lod1: are there phenoentries that equal LL and are LOD = 1
    N1_llim1 <-                   sum(phenofile$trait  < lower_limit)
    N1_llim2 <- if(N1_llim1 > 0L) sum(phenofile$trait  < lower_limit & phenofile$outsideLOD != 2L) else 0L
    N1_llim3 <-                   sum(phenofile$trait  > lower_limit & phenofile$outsideLOD == 2L)
    N1_llim4 <-                   sum(phenofile$trait == lower_limit & phenofile$outsideLOD == 0L)
    N1_llim_lod1 <-               sum(phenofile$trait == lower_limit & phenofile$outsideLOD == 1L)
    N1_llim_lod2 <-               sum(phenofile$trait == lower_limit & phenofile$outsideLOD == 2L)
  }
  
  
  ### checking upper limit
  
  N1_lod0 <- sum(phenofile$outsideLOD == 0L)
  if(is.na(upper_limit)){
    N1_ulim_lod1 <-0L
    N1_ulim_lod0 <- 0L
    N1_ulim4 <- 0L
  } else {
    # ulim1: are there pheno entries > UL
    # ulim2: are there phenos > UL not marked as LOD0
    # ulim3: are there lod0 entries that fall below UL
    # ulim4: are there phenoentries that equal UL, yet are LOD = 2
    N1_ulim1 <-                   sum(phenofile$trait  > upper_limit)
    N1_ulim2 <- if(N1_ulim1 > 0L) sum(phenofile$trait  > upper_limit & phenofile$outsideLOD != 0L) else 0L
    N1_ulim3 <-                   sum(phenofile$trait  < upper_limit & phenofile$outsideLOD == 0L)
    N1_ulim4 <-                   sum(phenofile$trait == upper_limit & phenofile$outsideLOD == 2L)
    N1_ulim_lod1 <-               sum(phenofile$trait == upper_limit & phenofile$outsideLOD == 1L)
    N1_ulim_lod0 <-               sum(phenofile$trait == upper_limit & phenofile$outsideLOD == 0L)
  }

  
  ### checking for reversal
  
  if(N1_lod2 > 0L) {
    # Condition 1: All values at LOD2 (lower limit) must be equal to or above values of LOD1
    reversal_low <- min(phenofile$trait[phenofile$outsideLOD == 2L]) >= max(phenofile$trait[phenofile$outsideLOD == 1L])
  }
  
  if(N1_lod0 > 0L) {
    # Condition 2: All values at LOD0 (upper limit) must be equal to or below values of LOD1
    reversal_hig <- max(phenofile$trait[phenofile$outsideLOD == 0L]) <= min(phenofile$trait[phenofile$outsideLOD == 1L])
    error_reversal <- if(N1_lod2 > 0L) reversal_low & reversal_hig else reversal_hig  
  } else {
    error_reversal <- if(N1_lod2 > 0L) reversal_low else FALSE
  }

  
  ### creating script output
  print("", quote = FALSE)
  print(" *** Results of lodQC ***", quote = FALSE)
  write.table(" *** Results of lodQC ***", paste0(filedirectory, "/", outputfile, ".log"), append = FALSE,
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  con_log <- file(description = paste0(filedirectory, "/", outputfile, ".log"), open = "a")
  
  out_sets <- data.frame(User = c("File", "phenotype", "lower_limit", "upper_limit"),
                         settings = c(outputfile, pheno_name, lower_limit, upper_limit),
                         stringsAsFactors = FALSE)
  print.data.frame(out_sets, right = TRUE, row.names = FALSE)
  print("", quote = FALSE)
  print("  *  The table below counts entries below, within and above LOD by", quote = FALSE)
  print("     using the 'outsideLOD' column in 'phenofile'.", quote = FALSE)
  print("     The lower_limit and upper_limit (see above) are defined by", quote = FALSE)
  print("     the user and used only to check if the 'outsideLOD' data", quote = FALSE)
  print("     matches the expected limits of detection", quote = FALSE)
  write.table("User settings", con_log, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  write.table(out_sets, con_log, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  write.table(c("",
                "  *  The table below counts entries below, within and above LOD by using the 'outsideLOD' column in 'phenofile'.",
                "     The lower_limit and upper_limit (see above) are defined by the user and used only to check if the 'outsideLOD' data matches the expected limits of detection."),
              con_log, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  
  
  ### creating summary statistic output table
  
  out_sumstat <- data.frame(Summary = c("Total", "   below LOD", "   within LOD", "   above LOD"),
                            empty = character(4),
                            N = c(N0_all, N0_lod2, N0_lod1, N0_lod0), valid = integer(4),
                            min = numeric(4), p5 = numeric(4), p25 = numeric(4),
                            p50 = numeric(4), p75 = numeric(4), p95 = numeric(4), max = numeric(4),
                            mean = numeric(4), sd = numeric(4), stringsAsFactors = FALSE)
  calc_sumstat <- function(LODval, DD, intempl) {
    useval <- if(is.na(LODval)) { rep(TRUE, nrow(DD)) } else { DD$outsideLOD == LODval }
    intempl$valid[1] <- sum(useval)
    
    if(intempl$valid[1] == 0L) {
      intempl[1, 5:13] <- NA
    } else {
      intempl[1, 5:11] <- quantile(DD[useval, 1],
                                  probs = c(0, 0.05, 0.025, 0.5, 0.75, 0.95, 1),
                                  names = FALSE)
      intempl$mean[1] <- mean(DD[useval, 1])
      intempl$sd[1]   <-   sd(DD[useval, 1])
      intempl[1, 5:13] <- round(intempl[1, 5:13], digits = 4)
    }
    return(intempl[1,])
  }
  
  out_sumstat[1, ] <- calc_sumstat(LODval = NA, DD = phenofile, intempl = out_sumstat[1, ])
  out_sumstat[2, ] <- calc_sumstat(LODval = 2L, DD = phenofile, intempl = out_sumstat[2, ])
  out_sumstat[3, ] <- calc_sumstat(LODval = 1L, DD = phenofile, intempl = out_sumstat[3, ])
  out_sumstat[4, ] <- calc_sumstat(LODval = 0L, DD = phenofile, intempl = out_sumstat[4, ])
  

  print("", quote = FALSE)
  print("", quote = FALSE)
  print.data.frame(out_sumstat[ , -2], right = FALSE, row.names = FALSE) # -2 removes empty column
  
  write.table(c("", ""), con_log, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  write.table(t(c("Summary", "", "N", "valid", "min", "5%", "25%", "median", "75%", "95%", "max", "mean", "standard deviation")),
              con_log, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  write.table(out_sumstat, con_log, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  
  
  ### adding notes if there are both censored and uncensored values at the upper and/or lower limits
  # we use the write error function for this, but do not generate warnings
  write_error <- function(error_message, generate_warning = reprint_warnings){
    if(generate_warning) warning(error_message, call. = FALSE)
    print(paste0("   ", error_message), quote = FALSE)
    return(paste0("   ", error_message))
  }
  
  if(N1_llim_lod2 + N1_ulim_lod0 > 0L){
    print("", quote = FALSE)
    print("", quote = FALSE)
    write.table(c("", ""), con_log, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)
    list_error <- NULL
    if(N1_llim_lod2 > 0L){
      list_error <- c(list_error,
                      write_error(paste0("There are ", N1_llim_lod2 + N1_llim_lod1 + N1_llim4,
                                         " entries whose value equals the 'lower_limit' (", lower_limit, ")."),
                                  generate_warning = FALSE),
                      write_error( paste("   Of these, 'outsideLOD' marks", N1_llim_lod1, "entries as within LOD; and",
                                         N1_llim_lod2, "as entries below LOD."),
                                  generate_warning = FALSE))
    }
    if(N1_ulim_lod0 > 0L){
      list_error <- c(list_error,
                      write_error(paste0("There are ", N1_ulim_lod1 + N1_ulim_lod0 + N1_ulim4,
                                         " entries whose value equals the 'upper_limit' (", upper_limit, ")."),
                                  generate_warning = FALSE),
                      write_error( paste("   Of these, 'outsideLOD' marks", N1_ulim_lod1, "entries as within LOD; and",
                                         N1_ulim_lod0, "entries as above LOD."),
                                  generate_warning = FALSE))
    }
    write.table(list_error, con_log, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)  
  }
  
  
  ### Outputting warnings
  print("", quote = FALSE)
  print("", quote = FALSE)  
  print(" * WARNINGS *", quote = FALSE)  
  write.table(c("", "", " * WARNINGS *"), con_log, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  error_found <- error_reversal # takes value from the reversal
  list_error <- NULL
  

  # error relating to missing phenotypes
  if(N1_lod2 == 0L){
    list_error <- c(list_error,
                    if(N0_lod2 == 0L) { write_error("outsideLOD marks 0 entries as below LOD")
                    } else { write_error("outsideLOD marks 0 valid phenotypes as below LOD") })
  }
  if(N1_lod0 == 0L){
    list_error <- c(list_error,
                    if(N0_lod0 == 0L) { write_error("outsideLOD marks 0 entries as above LOD")
                    } else { write_error("outsideLOD marks 0 valid phenotypes as above LOD") })
  }
  
  
  # if outsideLOD values are reversed, all other warnings and tests become redundant
  if(error_reversal){
    list_error <- c(list_error,
                    write_error("Column 'outsideLOD' uses inverted values. To run analysis correctly, use 2 for values below LOD and 0 for values above.",
                                generate_warning = FALSE)) # no warning necessary, at this will automatically terminate the function
  } else {
    # errors relating to the lower limit
    if(is.na(lower_limit)){
      if(N1_lod2  > 0L){ list_error <- c(list_error, write_error(paste("User specified no 'lower_limit', but 'outsideLOD' marks", N1_lod2, "entries as below LOD"))) }
    } else {
      if(N1_llim1 > 0L){
        list_error <- c(list_error, write_error(paste(N1_llim1, "phenotype entries have values below 'lower_limit'"),
                                                generate_warning = FALSE) )
      }
      if(N1_llim2 + N1_llim3 + N1_llim4 > 0L){
        error_found <- TRUE
        if(N1_llim2 > 0L){ list_error <- c(list_error, write_error(paste(N1_llim2, "phenotype entries are below   'lower_limit', but are *NOT* marked as such by 'outsideLOD'"))) }
        if(N1_llim3 > 0L){ list_error <- c(list_error, write_error(paste(N1_llim3, "phenotype entries *EXCEED*    'lower_limit', but are *NOT* marked as such by 'outsideLOD'"))) }
        if(N1_llim4 > 0L){ list_error <- c(list_error, write_error(paste(N1_llim4, "phenotype entries have a value equal to the 'lower_limit' but are marked as at the 'upper_limit' by 'outsideLOD'"))) }
      }
    }
    
    # error relating to the upper limit
    if(is.na(upper_limit)){
      if(N1_lod0  > 0L){ list_error <- c(list_error, write_error(paste("User specified no 'upper_limit', but 'outsideLOD' marks", N1_lod0, "entries as above LOD"))) }
    } else {
      if(N1_ulim1 > 0L){
        list_error <- c(list_error, write_error(paste(N1_ulim1, "phenotype entries have values above 'upper_limit'"),
                                                generate_warning = FALSE) )
      }
      if(N1_ulim2 + N1_ulim3 + N1_ulim4 > 0L){
        error_found <- TRUE
        if(N1_ulim2 > 0L){ list_error <- c(list_error, write_error(paste(N1_ulim2, "phenotype entries exceed      'upper_limit', but are *NOT* marked as such by 'outsideLOD'"))) }
        if(N1_ulim3 > 0L){ list_error <- c(list_error, write_error(paste(N1_ulim3, "phenotype entries are *BELOW* 'upper_limit', but are *NOT* marked as such by 'outsideLOD'"))) }
        if(N1_ulim4 > 0L){ list_error <- c(list_error, write_error(paste(N1_ulim4, "phenotype entries have a value equal to the 'upper_limit' but are marked as at the 'lower_limit' by 'outsideLOD'"))) }
      }
    }
    
    
    ### check for GAP between limit and extreme values
    # this must happen here, since the detection limit cannot be changed until above checks have run
    # also won't trigger if error_reversal is TRUE
    
    # checking lower limit
    temp_descr <- "Lower detection limit"
    if(is.na(lower_limit)) {
      if (N1_lod2 > 0L) {
        lower_limit <- max(phenofile$trait[phenofile$outsideLOD == 2L])
        temp_descr <- "Largest value for left-censored data"
      } }
    
    if(!is.na(lower_limit)) {
      # if no lower_limit and outsideLOD == 2L, this is skipped
      # lower limit is compared to min - (quantile[0.25] - min)
      # i.e. if it's the same distance from the minimum as the 0.25 quantile,
      # but in the other direction
      if (lower_limit <
            2 * min(phenofile$trait[phenofile$outsideLOD == 1L]) - quantile(phenofile$trait[phenofile$outsideLOD == 1L], 0.25)) {
        error_found <- TRUE
        list_error <- c(list_error, write_error(
          paste0(temp_descr, " (", signif(lower_limit, 4),
                 ") is much lower than smallest phenotype value within LOD (",
                 signif(min(phenofile$trait[phenofile$outsideLOD == 1L]), 4), ")")))
      } }
    
    
    # checking upper limit
    temp_descr <- "Upper detection limit"
    if(is.na(upper_limit)) {
      if (N1_lod0 > 0L) {
        upper_limit <- min(phenofile$trait[phenofile$outsideLOD == 0L])
        temp_descr <- "Smallest value for right-censored data"
      } }
    
    if(!is.na(upper_limit)) {
      # if no upper_limit and outsideLOD == 0L, this is skipped
      if (upper_limit >
            2 * max(phenofile$trait[phenofile$outsideLOD == 1L]) - quantile(phenofile$trait[phenofile$outsideLOD == 1L], 0.75)) {
        error_found <- TRUE
        list_error <- c(list_error, write_error(
          paste0(temp_descr, " (", signif(upper_limit, 4),
                 ") is much higher than largest phenotype value within LOD (",
                 signif(max(phenofile$trait[phenofile$outsideLOD == 1L]), 4), ")")))
      } }
  }
  
  if(is.null(list_error)){
    print(" - none", quote = FALSE)
    write.table(" - none", con_log, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)
  } else {
    write.table(list_error, con_log, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)  
  }
  
  close(con_log)  
  

  print("", quote = FALSE)
  print("", quote = FALSE)  
  if(error_reversal) stop("Column 'outsideLOD' uses inverted values. To run analysis correctly, use 2 for values below LOD and 0 for values above.")
  if(stop_if_error & error_found) stop("Column 'outsideLOD' does not match the specified limits")
  return(invisible(NULL))
}


lod_GWAS <- function(phenofile, pheno_name,
                    basic_model = NULL, dist = "gaussian",
                    mapfile, genofile, outputfile,
                    filedirectory = getwd(),
                    outputheader = "QCGWAS",
                    gzip_output = TRUE,
                    lower_limit = NA,
                    upper_limit = NA){
  geno_header <- TRUE # currently unused
  
  ### Checking input
  
  if(!is.character(filedirectory) | length(filedirectory) != 1L) stop("'filedirectory' must be a character-string containing one element")
  if(!file.exists(filedirectory)) stop("Cannot find 'filedirectory'")
  
  if(!is.character(outputfile) | length(outputfile) != 1L) stop("'outputfile' must be a character-string containing one element")
  if(!is.logical(gzip_output) | length(gzip_output) != 1L | is.na(gzip_output)) stop("Argument 'gzip_output' must be TRUE or FALSE")
  if(gzip_output){
    if(substr(outputfile, nchar(outputfile) - 2L, nchar(outputfile)) != ".gz"){
      outputfile <- paste0(outputfile, ".gz")
    }
  }
  #if(file.exists(paste0(filedirectory, "/", outputfile))) stop(paste0("Working directory already contains a file named '", outputfile, "'"))
  
  use_model <- !is.null(basic_model)
  if(use_model){
    if(!is.character(basic_model) | length(basic_model) != 1L) stop("'basic_model' must be a character-string containing one element")
  }
  
  # dist is not checked, as it can be a distribution as well
  if(is.list(dist)){
    if(!survregDtest(dlist = dist, verbose = FALSE)) stop("'dist' is not a valid survival distribution. Use function survregDtest to determine why.")
  } else {
    input_OK <- is.character(dist) & length(dist) == 1L
    if(input_OK) { input_OK <- dist %in% c("weibull", "exponential", "gaussian", "logistic","lognormal", "loglogistic") }
    if(!input_OK) stop("'dist' must either be a valid survival distribution, or a one of the following character strings: 'weibull', 'exponential', 'gaussian', 'logistic','lognormal', 'loglogistic'")
  }
  
  input_OK <- is.character(outputheader) & length(outputheader) == 1L
  if(input_OK) input_OK <- outputheader %in% c("QCGWAS", "oldQCGWAS", "GWAMA", "PLINK", "META", "GenABEL", "standard")
  if(!input_OK) stop("'outputheader' must be one of 'QCGWAS', 'oldQCGWAS', 'GWAMA', 'PLINK', 'META', 'GenABEL'")
  
  # Checking genofile (goes first because we don't have to load it)
  if(is.data.frame(genofile) | is.matrix(genofile)) stop("For technical reasons, the 'genofile' argument can only be a filename, not a dataset")
  if(!is.character(genofile) | length(genofile) != 1L) stop("'genofile' is not a single filename")
  if(!file.exists(paste0(filedirectory, "/", genofile))) stop(paste0("Cannot find genofile '", genofile, "'"))

  #if(!is.logical(geno_header) | length(geno_header) != 1L | any(is.na(geno_header))) stop("Argument 'geno_header' must be a single boolean")  
  #if(!geno_header) print("User specified that genotype file does not sample IDs. The function will assume that the values in the phenotype data already match the order of samples in the genotype file.", quote = FALSE)
  
  
  # Checking the phenotype file
  loadandreturn <- function(filename, filedirectory, header = FALSE){
    
    if(is.data.frame(filename) | is.matrix(filename)) return(filename)
    
    if(!is.character(filename) | length(filename) != 1L) stop("One of the filename arguments is neither a dataset nor a filename")
    
    if(!file.exists(paste0(filedirectory, "/", filename))) stop(paste0("Cannot find file '", filename, "' in the working directory"))
    
    return(read.table(file = paste0(filedirectory, "/", filename), stringsAsFactors = FALSE, header = header, comment.char = ""))
  }
  
  if(!is.character(pheno_name) | length(pheno_name) != 1L) stop("'pheno_name' must be a character-string containing one element")
  
  phenofile <- loadandreturn(filename = phenofile, filedirectory = filedirectory, header = TRUE)
  pheno_colnames <- colnames(phenofile)  
  if(any(duplicated(pheno_colnames))) stop("Phenotype file contains duplicated column names")
  if(!all(c("FID", "IID", "outsideLOD") %in% pheno_colnames)) stop(paste("Phenotype file does not contain FID, IID and/or outsideLOD columns. Columns present in phenotype file are:", paste(if(length(pheno_colnames) < 20L) pheno_colnames else pheno_colnames[1:20], collapse = ", ")))
  
  pheno_cols <- which(pheno_colnames == pheno_name)
  if(length(pheno_cols) == 1L) {
    if(any(pheno_colnames[-pheno_cols] == "trait")) stop(paste0("Cannot assign column '", pheno_name, "' as trait because there already is a column named 'trait'"))
  #  colnames(phenofile)[pheno_cols] <- "trait"
  } else{
    if(length(pheno_cols) == 0L) { stop(paste0("Dataset 'phenofile' contains no column named '", pheno_name, "'. Columns present in phenotype file are: ", paste(if(length(pheno_colnames) < 20L) pheno_colnames else pheno_colnames[1:20], collapse = ", ")))
    } else { stop(paste0("Dataset 'phenofile' contains multiple columns named '", pheno_name, "'")) }
  }

  # NA values in outsideLOD are ignored if the correspond trait value is also NA
  input_OK <- is.integer(phenofile$outsideLOD) & !any(is.na(phenofile$outsideLOD) & !is.na(phenofile[ , pheno_cols]))
  if(input_OK) { input_OK <- all(phenofile$outsideLOD > -1L, na.rm = TRUE) & all(phenofile$outsideLOD < 3L, na.rm = TRUE) }
  if(!input_OK) stop("The column 'outsideLOD' can only contain integer values 0, 1 and 2")
  
  # checking if basic model factors are present
  if(use_model){
    model_elements <- unlist(strsplit(basic_model, split = "[:+:]"))
    if(any(model_elements %in% c("trait", pheno_name, "outsideLOD", "SNP"))) stop("The argument 'basic_model' can only contain covariates. Do not include 'trait', 'outsideLOD', 'SNP' or the name of the phenotype.")
    if(any(!model_elements %in% pheno_colnames)) { stop(paste("Covariates", paste(model_elements[!model_elements %in% pheno_colnames], collapse = ", "), "not present in dataset")) }
  } else model_elements <- NULL
  rm(pheno_colnames)

  
  # checking LOD
  input_OK <- is.vector(lower_limit) & length(lower_limit) == 1L
  if(input_OK) { input_OK <- is.numeric(lower_limit) | is.na(lower_limit) }
  if(!input_OK) stop("'lower_limit' must be either a numeric or missing value")
  
  input_OK <- is.vector(upper_limit) & length(upper_limit) == 1L
  if(input_OK) { input_OK <- is.numeric(upper_limit) | is.na(upper_limit) }
  if(!input_OK) stop("'upper_limit' must be either a numeric or missing value")

  lod_QC(phenofile = phenofile, pheno_name = pheno_name,
        lower_limit = lower_limit, upper_limit = upper_limit,
        stop_if_error = FALSE, reprint_warnings = TRUE,
        filedirectory = filedirectory, outputfile = outputfile)
  # phenotype column renamed after the lod_QC
  colnames(phenofile)[pheno_cols] <- "trait"

  # Loading & checking map file
  mapfile <- loadandreturn(filename =  mapfile, filedirectory = filedirectory)
  if(ncol(mapfile) != 4L) stop("'mapfile' must have 4 columns")
  colnames(mapfile) <- c("chr","SNP","cm","pos")  
  

  ### Start of analysis
  
  #con_geno <- file(paste0(filedirectory, "/", genofile), open = "rt") # this doesn't seem to work
  print("", quote = FALSE)
  con_geno <- paste0(filedirectory, "/", genofile)
  if(geno_header){
    dosage_header <- read.table(con_geno, header=FALSE, nrows=1L, stringsAsFactors=FALSE)
    dosage_ID <- data.frame(FID = t(dosage_header[ , 2+2*1:((length(dosage_header) - 3) / 2 )]),
                            IID = t(dosage_header[ , 3+2*1:((length(dosage_header) - 3) / 2 )]),
                            stringsAsFactors = FALSE)
    rm(dosage_header)
    
    N_samples <- sum(paste(phenofile$FID, phenofile$IID, sep = "_") %in% paste(dosage_ID$FID, dosage_ID$IID, sep = "_"))
    N_samples_missing_geno  <- nrow(phenofile) - N_samples
    N_samples_missing_pheno <- nrow(dosage_ID) - N_samples
    
    if(N_samples < 5L) stop("Less than 5 matching IDs found between geno- and phenotype files")
    if(N_samples_missing_geno  > 0L) print(paste(N_samples_missing_geno , "phenotyped samples genotype data"),
                                           quote = FALSE)
    if(N_samples_missing_pheno > 0L) print(paste(N_samples_missing_pheno, "genotyped samples do not appear in the phenotype data"),
                                           quote = FALSE)
    nrlines <- length(readLines(con_geno)) - 1L
  } else { # currently not used - code incomplete
    #dosage_ID <- data.frame(FID = phenofile$FID,
    #                        IID = phenofile$IID,
    #                        stringsAsFactors = FALSE)
    #nrlines <- length(readLines(con_geno))
    stop("This code should never run - geno_header must be TRUE")
  }
  
  # a wrapper function run survivalf with sapply
  apply_survivalf <- function(iter, dosage, ...) return(survivalf(SNP = dosage[ , iter+2], ...))
  N_SNPs_not_in_mapfile <- 0L
  nrrows <- 5000L
  k <- 0L
  out_results <- data.frame(MARKER = character(nrlines),
                            CHR = if(is.integer(mapfile$chr)) integer(nrlines) else character(nrlines),
                            POSITION = integer(nrlines),
                            OTHER_ALL = character(nrlines), EFFECT_ALL = character(nrlines),
                            N_TOTAL = integer(nrlines), N_VALID = integer(nrlines),
                            EFF_ALL_FREQ = numeric(nrlines),
                            EFFECT = numeric(nrlines), STDERR = numeric(nrlines),
                            PVALUE = numeric(nrlines), IMP_QUALITY= numeric(nrlines),
                            stringsAsFactors = FALSE)
  print("", quote = FALSE)
  while (k < nrlines) {
    if(k!=0L) {
      print(paste0("   reading SNPs ", k+1L, " to ", if(nrlines > k+nrrows) k+nrrows else nrlines), quote = FALSE)
      flush.console()
    }
    dosage <- read.table(con_geno, header=FALSE, skip=if(geno_header) 1L+k else k, nrows=nrrows, stringsAsFactors=FALSE)
    SNPdetails <- data.frame(dosage[,1:3], order=c(1:nrow(dosage)))
    names(SNPdetails) <- c("SNP", "A1", "A2", "order")
    dosage <- dosage[ , 4:ncol(dosage)]
    
    # determining dosage type
    if(k==0L){
      if(!is.numeric(dosage[,1])) stop(
        if(geno_header) { "Non-numeric data in the genotype dosages"
                          } else { "Non-numeric data in the genotype dosages - user specified no header" } )
      if (ncol(dosage) == nrow(dosage_ID)) {
        print("Dosage data format is 1: dosages", quote=FALSE)
        dtype <- 1L
      } else if (ncol(dosage) == 2*nrow(dosage_ID)) {
        print("Dosage data format is 2: two probabilities", quote=FALSE)
        dtype <- 2L
      } else if (ncol(dosage) == 3*nrow(dosage_ID)) {
        print("Dosage data format is 3: three probabilities",quote=FALSE)
        dtype <- 3L
      } else {
        #close(con_geno)
        stop(if(geno_header) "Dosage data format is not recognized" else "The number of genotyped samples does not match the number of samples in phenofile")
      }
      print(paste0("   reading SNPs 1 to ", if(nrlines > nrrows) nrrows else nrlines), quote = FALSE)
      flush.console() # this message is in fact false, since the dosages have already been read. But otherwise the console output looks odd
    }
    
    if(dtype == 1L) {
      dosage <- data.frame(t(dosage)) }
    if(dtype == 2L) {
      dosage <- data.frame(t(dosage[ , 2*1:(ncol(dosage)/2)-1]*2 + dosage[,2*1:(ncol(dosage)/2)])) }
    if(dtype == 3L) {
      dosage <- data.frame(t(dosage[ , 3*1:(ncol(dosage)/3)-2]*2 + dosage[ , 3*1:(ncol(dosage)/3)-1])) }
    
    
    names(dosage) <- SNPdetails[ , 1]
    dosage <- merge(data.frame(dosage_ID, dosage),
                    phenofile[ , c("FID", "IID", "trait", "outsideLOD", model_elements)],
                    by = c("FID", "IID"), all.x = TRUE, all.y = FALSE)
    if(k==0L){
      phenotemp <- dosage[ , c("trait", "outsideLOD", model_elements)]
      
      list_missing <- is.na(phenotemp$trait)
      if((length(list_missing) - sum(list_missing)) < 5L) {
        stop("There are insufficient non-missing phenotype values for the genotyped samples") }
      if(use_model){
        for(lmi in 1:length(model_elements)){
          list_missing <- list_missing | is.na(phenotemp[ , model_elements[lmi]]) }
        if((length(list_missing) - sum(list_missing)) < 5L) {
          stop("There are insufficient non-missing covariates for the genotyped samples") }
      }
      rm(list_missing)
    }

    N_SNPs_not_in_mapfile <- N_SNPs_not_in_mapfile + sum(!SNPdetails$SNP %in% mapfile$SNP)
    SNPdetails <- merge(SNPdetails, mapfile, by = "SNP", all.x = TRUE, all.y = FALSE)
    SNPdetails <- SNPdetails[order(SNPdetails[,"order"]), ]
    
    results_SNPs <- t(sapply(X = 1:nrow(SNPdetails), FUN = apply_survivalf,
                             dosage = dosage, data = phenotemp, basic_model = basic_model, dist = dist))
    out_results[(k+1):(k+nrow(SNPdetails)), ] <- cbind(SNPdetails[ , c("SNP", "chr", "pos", "A2", "A1")], results_SNPs)
    k <- k + nrrows
  }

  ### adding file header
  
  ### generating file header
  if(outputheader == "QCGWAS"){
    colnames(out_results) <- c("MARKER", "CHR", "POSITION",
                           "OTHER_ALL", "EFFECT_ALL",
                           "N_TOTAL", "N_VALID", "EFF_ALL_FREQ",
                           "EFFECT", "STDERR", "PVALUE", "IMP_QUALITY")
  } else if(outputheader == "oldQCGWAS") {
    colnames(out_results) <- c("MARKER", "CHR", "POSITION",
                           "ALLELE2", "ALLELE1",
                           "N_TOTAL", "N_VALID", "FREQLABEL",
                           "EFFECT", "STDERR", "PVALUE", "IMP_QUALITY")
  } else if(outputheader == "GWAMA"){
    colnames(out_results) <- c("MARKER", "CHR", "POSITION",
                           "NEA", "EA",
                           "N", "N_valid", "EAF",
                           "BETA", "SE", "P", "IMP_QUALITY")
  } else if(outputheader == "PLINK") {    
    colnames(out_results) <- c("SNP", "CHR", "BP",
                           "A2", "A1",
                           "N", "N_valid", "EFF_ALL_FREQ",
                           "BETA", "SE", "P", "IMP_QUALITY")
  } else if(outputheader == "META") {
    colnames(out_results) <- c( "rsid",   "chr", "pos",
                            "allele_A", "allele_B",
                            "N", "N_valid", "EFF_ALL_FREQ",
                            "beta", "se", "P_value", "info") 
  } else if(outputheader == "GenABEL") {
    colnames(out_results) <- c("name", "chromosome", "position",
                           "allele2", "allele1",
                           "n", "n_valid", "effallelefreq",
                           "beta", "sebeta", "p", "impquality")
  } else { # default option
    colnames(out_results) <- c("Markerid", "Chr", "Position",
                           "Other_allele", "Coded_allele",
                           "N_total", "N_valid", "Coded_allele_freq",
                           "Beta", "SE", "P", "Imp_info")
  }
  
  write.table(out_results,
              if(gzip_output) gzfile(paste0(filedirectory, "/", outputfile)) else paste0(filedirectory, "/", outputfile),
              col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
  

  print("", quote = FALSE)
  if(N_SNPs_not_in_mapfile > 0L) print(paste(N_SNPs_not_in_mapfile, "SNPs not found in mapfile"), quote = FALSE)
  print("", quote = FALSE)
  return(invisible(NULL))
}
