summaryData.setup <- function(summary.files, pathway, family, reference, lambda, 
                              ncases, ncontrols, nsamples, options){
  
  start.time <- date()
  
  # validate the format of main inputs
  validate.summary.input(summary.files, pathway, family, reference, lambda, 
                         ncases, ncontrols, nsamples)
  
  # reformat reference to convert factor to character
  reference <- reformat.reference.path(reference)
  
  # merge and reset options
  options <- options.setup(options, family, lambda, ncases, ncontrols, nsamples)
  
  # load definition of pathway
  pathway <- load.pathway.definition(pathway, options)
  
  # load and check summary statistics
  sum.stat <- load.summary.statistics(summary.files, family, pathway$SNP, options)
  
  # deleted snps and their reason
  deleted.snps <- data.frame(SNP = NULL, reason = NULL, comment = NULL, stringsAsFactors = FALSE)
  deleted.genes <- data.frame(Gene = NULL, reason = NULL, stringsAsFactors = FALSE)
  exc.snps <- intersect(pathway$SNP, options$excluded.snps)
  exc.snps <- setdiff(exc.snps, options$selected.snps)
  deleted.snps <- update.deleted.snps(deleted.snps, exc.snps, reason = "RM_BY_SNP_NAMES", comment = "")
  pathway <- update.pathway.definition(pathway, exc.snps)
  sum.stat <- update.sum.stat(sum.stat, exc.snps)
  
  # delete SNPs in specific regions
  exc.reg <- find.snps.in.regions(sum.stat$stat, options)
  exc.snps <- exc.reg$exc.snps
  deleted.snps <- update.deleted.snps(deleted.snps, exc.snps, reason = "RM_BY_REGIONS", comment = exc.reg$comment)
  pathway <- update.pathway.definition(pathway, exc.snps)
  sum.stat <- update.sum.stat(sum.stat, exc.snps)
  
  # update with valid/available SNPs
  exc.snps <- setdiff(pathway$SNP, sum.stat$snps.in.study)
  deleted.snps <- update.deleted.snps(deleted.snps, exc.snps, reason = "NO_SUM_STAT", comment = "")
  pathway <- update.pathway.definition(pathway, exc.snps)
  sum.stat <- update.sum.stat(sum.stat, exc.snps)
  
  # load SNPs and their reference and effect alleles in reference genotype
  allele.info <- load.reference.allele(reference, pathway, options)
  ref.snps <- allele.info$SNP
  
  # update with valid/available SNPs
  exc.snps <- setdiff(pathway$SNP, ref.snps)
  deleted.snps <- update.deleted.snps(deleted.snps, exc.snps, reason = "NO_REF", comment = "")
  pathway <- update.pathway.definition(pathway, exc.snps)
  sum.stat <- update.sum.stat(sum.stat, exc.snps)
  allele.info <- update.allele.info(allele.info, exc.snps)
  ref.snps <- update.ref.snps(ref.snps, exc.snps)
  
  # filter out SNPs with conflictive allele information in summary statistics and reference
  exc.snps <- filter.conflictive.snps(sum.stat, allele.info, options)
  
  # update with valid/available SNPs
  deleted.snps <- update.deleted.snps(deleted.snps, exc.snps, reason = "CONF_ALLELE_INFO", comment = "")
  pathway <- update.pathway.definition(pathway, exc.snps)
  sum.stat <- update.sum.stat(sum.stat, exc.snps)
  allele.info <- update.allele.info(allele.info, exc.snps)
  ref.snps <- update.ref.snps(ref.snps, exc.snps)
  
  # load genotypes in reference
  ref.geno <- load.reference.geno(reference, pathway$SNP, options)
  exc.snps <- setdiff(pathway$SNP, colnames(ref.geno))
  deleted.snps <- update.deleted.snps(deleted.snps, exc.snps, reason = "SNP_CONST", comment = "")
  pathway <- update.pathway.definition(pathway, exc.snps)
  sum.stat <- update.sum.stat(sum.stat, exc.snps)
  allele.info <- update.allele.info(allele.info, exc.snps)
  ref.snps <- update.ref.snps(ref.snps, exc.snps)
  
  ar <- align.reference(ref.geno, allele.info, options)
  ref.geno <- ar$ref.geno
  allele.info <- ar$allele.info
  
  # SNP filtering based on options
  filtered.data <- filter.reference.geno(ref.geno, pathway, sum.stat, options)
  filtered.markers <- filtered.data$deleted.snps
  filtered.genes <- filtered.data$deleted.genes
  
  # update with valid/available SNPs
  exc.snps <- filtered.markers$SNP
  exc.genes <- filtered.genes$Gene
  deleted.snps <- update.deleted.snps(deleted.snps, exc.snps, 
                                      reason = filtered.markers$reason, 
                                      comment = filtered.markers$comment)
  deleted.genes <- update.deleted.genes(deleted.genes, exc.genes, filtered.genes$reason)
  pathway <- update.pathway.definition(pathway, exc.snps, exc.genes)
  sum.stat <- update.sum.stat(sum.stat, exc.snps)
  allele.info <- update.allele.info(allele.info, exc.snps)
  ref.snps <- update.ref.snps(ref.snps, exc.snps)
  ref.geno <- update.ref.geno(ref.geno, exc.snps)
  
  # estimate P and SE if they are not provided by users
  sum.stat <- complete.sum.stat(sum.stat, options)
  
  # redefine the pathway by introducing independent group to reduce computational burden
  pathway <- split.pathway(pathway, allele.info, options$group.gap)
  
  # recover the summary statistics
  norm.stat <- recover.stat(sum.stat, pathway, ref.geno, allele.info, options)
  
  if(!options$keep.geno){
    rm(ref.geno)
    gc()
    ref.geno <- NULL
  }
  
  # trim the information of deleted SNPs
  deleted.snps <- trim.deleted.snps(deleted.snps, options)
  
  msg <- paste0("Setup completed: ", date())
  if(options$print) message(msg)
  
  end.time <- date()
  setup.timing <- as.integer(difftime(strptime(end.time, "%c"), strptime(start.time, "%c"), units = "secs"))
  
  if(family == 'gaussian'){
    options$ncases <- NULL
    options$ncontrols <- NULL
  }
  
  if(options$meta){
    meta.stat <- meta(summary.files, lambda, sel.snps = unique(pathway$SNP), only.meta = options$only.meta)$meta.stat
  }else{
    meta.stat <- NULL
  }
  
  setup <- list(deleted.snps = deleted.snps, deleted.genes = deleted.genes, 
                options = options, meta.stat = meta.stat, allele.info = allele.info, 
                pathway = pathway, norm.stat = norm.stat, 
                ref.geno = ref.geno, setup.timing = setup.timing)
  
  if(options$save.setup){
    save(setup, file = options$path.setup)
    msg <- paste0("setup file has been saved at ", options$path.setup)
    message(msg)
  }
  
  setup
  
}
