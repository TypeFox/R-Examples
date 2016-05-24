
rawData.dataframe.setup <- function(formula, data, pathway, family, lambda, subset = NULL, options = NULL){
  
  start.time <- date()
  
  validate.family(family)
  validate.lambda.rawData(lambda)
  
  # merge and reset options
  options <- options.setup(options, family, lambda, NULL, NULL, NULL)
  
  # subset of data
  data <- data.subset(data, subset, options)
  
  # load definition of pathway
  pathway <- load.pathway.definition(pathway, options)
  
  # deleted snps and their reason
  deleted.snps <- data.frame(SNP = NULL, reason = NULL, comment = NULL, stringsAsFactors = FALSE)
  deleted.genes <- data.frame(Gene = NULL, reason = NULL, comment = NULL, stringsAsFactors = FALSE)
  exc.snps <- intersect(pathway$SNP, options$excluded.snps)
  exc.snps <- setdiff(exc.snps, options$selected.snps)
  deleted.snps <- update.deleted.snps(deleted.snps, exc.snps, reason = "RM_BY_SNP_NAMES", comment = "")
  pathway <- update.pathway.definition(pathway, exc.snps)
  
  # 
  exc.snps <- setdiff(pathway$SNP, colnames(data))
  deleted.snps <- update.deleted.snps(deleted.snps, exc.snps, reason = "NO_RAW_GENO", comment = "")
  pathway <- update.pathway.definition(pathway, exc.snps)
  
  #
  snps <- intersect(colnames(data), pathway$SNP)
  raw.geno <- data[, snps, drop = FALSE]
  null <- data[, setdiff(colnames(data), snps), drop = FALSE]
  
  ini <- data.parse(formula, null)
  rm(null, data)
  gc()
  null <- ini$null
  resp.var <- ini$resp.var
  comp.id <- which(complete.cases(null))
  if(length(comp.id) < nrow(null)){
    
    msg <- paste0(nrow(null) - length(comp.id), " samples are excluded due to missing covariates. ", length(comp.id), " samples are used")
    message(msg)
    
    null <- null[comp.id, ]
    
    #raw.geno <- raw.geno[comp.id, ]
    # the code above will double the memory consumption
    # use the codes below instead
    # the idea comes from @vc273 at 
    # http://stackoverflow.com/questions/10790204/how-to-delete-a-row-by-reference-in-r-data-table
    
    snps <- colnames(raw.geno)
    setDT(raw.geno)
    tmp <- data.table(V1 = raw.geno[[snps[1]]][comp.id])
    setnames(tmp, names(tmp), snps[1])
    for(s in snps[-1]){
      tmp[, s := raw.geno[[s]][comp.id], with = FALSE]
      raw.geno[, s := NULL, with = FALSE]
    }
    gc()
    raw.geno <- tmp
    class(raw.geno) <- "data.frame"
    rm(tmp)
    gc()
  }
  
  null <- validate.covar(null, resp.var)
  
  control.id <- which(null[, resp.var] == 0)
  
  msg <- paste("Filtering SNPs:", date())
  if(options$print) message(msg)
  
  uni.gene <- unique(pathway$Gene)
  for(g in uni.gene){
    gid <- which(pathway$Gene == g)
    
    if(length(gid) == 0){
      next
    }
    
    sid <- which(colnames(raw.geno) %in% pathway$SNP[gid])
    # SNP filtering based on options
    filtered.data <- filter.raw.geno(raw.geno[, sid, drop = FALSE], pathway[gid, , drop = FALSE], options, control.id)
    filtered.markers <- filtered.data$deleted.snps
    gc()
    
    # update with valid/available SNPs
    exc.snps <- filtered.markers$SNP
    deleted.snps <- update.deleted.snps(deleted.snps, exc.snps, 
                                        reason = filtered.markers$reason, 
                                        comment = filtered.markers$comment)
    pathway <- update.pathway.definition(pathway, exc.snps)
    raw.geno <- update.raw.geno(raw.geno, exc.snps)
    gc()
  }
  
  uni.chr <- unique(pathway$Chr)
  for(c in uni.chr){
    cid <- which(pathway$Chr == c)
    
    if(length(cid) == 0){
      next
    }
    
    sid <- which(colnames(raw.geno) %in% pathway$SNP[cid])
    # SNP filtering based on options
    filtered.data <- filter.raw.geno(raw.geno[, sid, drop = FALSE], pathway[cid, , drop = FALSE], options, control.id)
    filtered.markers <- filtered.data$deleted.snps
    filtered.genes <- filtered.data$deleted.genes
    gc()
    
    # update with valid/available SNPs
    exc.snps <- filtered.markers$SNP
    exc.genes <- filtered.genes$Gene
    deleted.snps <- update.deleted.snps(deleted.snps, exc.snps, 
                                        reason = filtered.markers$reason, 
                                        comment = filtered.markers$comment)
    deleted.genes <- update.deleted.genes(deleted.genes, exc.genes, filtered.genes$reason)
    pathway <- update.pathway.definition(pathway, exc.snps, exc.genes)
    raw.geno <- update.raw.geno(raw.geno, exc.snps)
    gc()
  }
  
  # calculate normal covariance and mean
  norm.stat <- generate.normal.statistics(resp.var, null, raw.geno, pathway, family, lambda)
  
  if(!options$keep.geno){
    rm(raw.geno)
    rm(null)
    gc()
    raw.geno <- NULL
    null <- NULL
  }else{
    class(raw.geno) <- 'data.frame'
  }
  
  # trim the information of deleted SNPs
  deleted.snps <- trim.deleted.snps(deleted.snps, options)
  
  msg <- paste0("Setup completed: ", date())
  if(options$print) message(msg)
  
  end.time <- date()
  setup.timing <- as.integer(difftime(strptime(end.time, "%c"), strptime(start.time, "%c"), units = "secs"))
  
  yx <- create.yx(resp.var, null)
  formula <- create.formula(resp.var, yx)
  
  setup <- list(deleted.snps = deleted.snps, deleted.genes = deleted.genes, 
                options = options, pathway = pathway, norm.stat = norm.stat, 
                formula = formula, yx = yx, raw.geno = raw.geno, 
                setup.timing = setup.timing)
  
  if(options$save.setup){
    save(setup, file = options$path.setup)
    msg <- paste0("setup file has been saved at ", options$path.setup)
    message(msg)
  }
  
  setup
  
}
