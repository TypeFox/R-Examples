alleleEnrichment <- function(
                      subsample.size,
                      preferred.alleles,
                      map, 
                      ped, 
                      sampleid.columnindex = 2,
                      mode = "max",
                      pheno.file = NULL,
                      na.samples = na.omit, 
                      na.snps = na.pass,
                      write.resultfiles = TRUE, 
                      ...
                    ) {

  args <- list(...)
  if(length(args) > 0 && !all((names(args) %in% c("num.perms", "search.seconds", "weights"))))
    stop("Invalid arguments supplied (in ...) or misspelled argument.\n")
  
  if(missing(subsample.size) || !is.numeric(subsample.size) || subsample.size <= 1)
    stop("The argument 'subsample.size' has to be a positive integer.\n")
    
  if(missing(preferred.alleles) || is.null(preferred.alleles) || (!is.list(preferred.alleles) && !is.atomic(preferred.alleles)) || is.null(names(preferred.alleles))) 
    stop("Argument 'preferred.alleles' has to be set and must be a named list or vector (names = SNP identifier, value = allele).\n")
  
  if(is.null(mode))
    stop("Argument 'mode' has to be set.\n")

  alleles <- alleleEnrichment.data.alleles(
               refalleles = preferred.alleles,
               map, 
               ped,
               pheno.file = pheno.file, 
               sampleid.columnindex, 
               na.samples = na.samples,
               na.snps = na.snps
             )

  if(ncol(alleles) <= subsample.size)
    stop("After NA removal, the remaining number of available samples is less than subsample.size\n")
  
  if(is.null(args$num.perms))
    num.perms <- 100 * length(preferred.alleles) * subsample.size
  else
    num.perms <- args$num.perms
  
  if(is.null(args$weights)) {
    weights <- as.list(rep(1, nrow(alleles)))
    names(weights) <- rownames(alleles)
  } else {
    weights <- args$weights
    weights <- unlist(weights)
    if(!is.atomic(weights) || length(weights) != nrow(alleles) || !all(names(weights) %in% rownames(alleles)) || !is.numeric(weights) || any(weights <= 0))
      stop("Invalid specification of weights. Must be a named list or vector according to enriched variables with numeric values > 0.\n")
    # map to interval [0,1] to prevent too large stepsize in the iterative simplex (leads to inaccurate results)
    weights <- as.list(weights / max(weights))
  }
  
  if(is.null(args$search.seconds))
    search.seconds <- 20
  else
    search.seconds <- args$search.seconds
  
  
  # make mode a list named by snps/phenotypes
  if(length(mode) == 1) {
    if(!is.character(mode) || (mode != "max" && mode != "balanced"))
      stop("Argument 'mode' has to be either 'max' or 'balanced' when specified as single element.\n")
    # set individual mode for each SNP + pheno by recycling mode argument
    mode <- lapply(rownames(alleles), function(name) unlist(mode))
    names(mode) <- rownames(alleles)
  } else {
    if(!is.list(mode) || is.null(names(mode)))
      stop("Argument 'mode' has to be of type 'list' when specified with multiple elements, and has to be named with SNPs / phentoypes.\n")
    if(length(mode) != nrow(alleles))
      stop(paste("Argument 'mode' has", length(mode), "elements, but needs", nrow(alleles), "elements.\n"))
    if(!all(names(mode) %in% rownames(alleles)))
      stop("Names of argument 'mode' do not match SNP and pehnotypes identifiers.\n")
    if(length(unique(lapply(mode, class))) != 1)
      stop("Elements of argument 'mode' have to be of the same type (e.g. cannot mix 'function' and 'character').\n")
  }

  # score having large amounts of risk alleles evenly distributed over SNPs (low standard deviation)
  # distributes evenly at its minimum
  deviation.weighted.minimize <- function(var.min, obs.varsum, weight)
    return(sqrt(sum(((var.min - obs.varsum)^2 * (1-weight)), na.rm = TRUE) / length(var.min)))

  # larger deviations do not contribute more than smaller. can be used for an unbiased maximization of SD
  deviation.weighted.maximize <- function(var.min, obs.varsum, weight)
    return(sum((abs(var.min - obs.varsum) * (1-weight)), na.rm = TRUE) / length(var.min))


  # allele count (sum) and half allele frequency (mean) in entire sample - also makes sense for most phentoypes
  var.sum <- apply(alleles, 1, sum)
  # required minimum sum in subsample
  var.min <- floor(var.sum * subsample.size / ncol(alleles))
  
  cat(paste("Choosing from", ncol(alleles), "genotyped samples (after NA exclusion)\n"))

  
  ######## optimize ########

  # use random trials only when all constraints are balanced    
  if(!all(mode == "balanced")) {
    
    sol <- alleleEnrichment.lp(
             subsample.size = subsample.size, 
             mode = mode, 
             alleles = alleles, 
             var.min = var.min, 
             weights = weights, 
             search.seconds = search.seconds
           )
    sol.varsum <- alleles %*% sol
    # quality is the summed differences with var.min
    sol.qual <- deviation.weighted.maximize(var.min, sol.varsum, 0)
    
  } else {
  
    sol <- alleleEnrichment.perm(
      subsample.size = subsample.size, 
      select.fun = deviation.weighted.minimize, 
      alleles = alleles, 
      var.min = var.min, 
      num.perms = num.perms
    )
    sol.varsum <- alleles %*% sol
    # quality is the summed differences with var.min
    sol.qual <- deviation.weighted.minimize(var.min, sol.varsum, 0)
  
  }
  
  
  ######## write summary files ########
  
  if(write.resultfiles) {
    writeLines(
      c(
        paste("Selected ", sum(sol), "out of", ncol(alleles), "samples with the following IDs:"),
        "",
        paste(colnames(alleles)[sol == 1], collapse = " "), 
        "",
        "Risk allele counts in the entire cohort:", 
        paste(rownames(alleles), var.sum, sep = "\t\t"),
        "",
        "Risk allele counts in the subsample:",
        paste(rownames(alleles), as.vector(sol.varsum), sep = "\t\t"),
        "",
        "Enrichment in percent of the original frequency:", 
        paste(rownames(alleles), round(as.vector(sol.varsum) / var.min * 100, digits = 2), sep = "\t\t")
      ),
      "alleleEnrichment_summary.txt"
    )
    writeLines(colnames(alleles)[sol == 1], "alleleEnrichment_sampleSelection.txt")
    write.table(
      alleles[, as.logical(sol)], 
      "alleleEnrichment_genotypes.txt", 
      quote = F, 
      sep = "\t"
    )
  }
  
  return(list(
      sample.selection = colnames(alleles)[sol == 1], 
      allele.counts = data.frame(variable = rownames(alleles), count = as.vector(sol.varsum)),
      enriched.percent = data.frame(variable = rownames(alleles), frequency = as.vector(sol.varsum) / var.min * 100),
      genotypes = alleles[, as.logical(sol)]
    ))

}
    
    