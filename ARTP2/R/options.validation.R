
options.validation <- function(options){
  
  if(is.null(options$id.str)){
    warning("options$id.str is NULL")
  }
  
  if(!any(options$method %in% 1:3)){
    stop("method should be 1 (AdaJoint), 2 (AdaJoint2), or 3 (ARTP)")
  }
  
  if(!is.numeric(options$huge.gene.size) || options$huge.gene.size < 0){
    stop("huge.gene.size should be non-negative integer")
  }
  
  if(!is.numeric(options$huge.chr.size) || options$huge.chr.size < 0){
    stop("huge.chr.size should be non-negative integer")
  }
  
  if(is.null(options$nperm)){
    stop("nperm cannot be NULL")
  }
  
  if(!is.null(options$nperm) && options$nperm < 1000){
    warning("nperm is too small")
  }
  
  if(!is.null(options$excluded.regions)){
    excluded.regions <- options$excluded.regions
    tmp <- (c("data.frame", "matrix") %in% class(excluded.regions))
    if(!any(tmp)){
      msg <- "options$excluded.regions should be either a data frame or a matrix"
      stop(msg)
    }else{
      if("matrix" %in% class(excluded.regions)){
        excluded.regions <- as.data.frame(excluded.regions)
      }
    }
    
    header1 <- c('Chr', 'Pos', 'Radius')
    header2 <- c('Chr', 'Start', 'End')
    if(all(header1 %in% colnames(excluded.regions))){
      excluded.regions <- excluded.regions[, header1, drop = FALSE]
      if(!any(c('integer', 'numeric') %in% class(excluded.regions$Pos))){
        msg <- 'options$excluded.regions should have integer Pos'
        stop(msg)
      }
      
      if(!any(c('integer', 'numeric') %in% class(excluded.regions$Radius))){
        msg <- 'options$excluded.regions should have integer Radius'
        stop(msg)
      }
    }else{
      if(all(header2 %in% colnames(excluded.regions))){
        excluded.regions <- excluded.regions[, header2, drop = FALSE]
        if(!any(c('integer', 'numeric') %in% class(excluded.regions$Start))){
          msg <- 'options$excluded.regions should have integer Start'
          stop(msg)
        }
        
        if(!any(c('integer', 'numeric') %in% class(excluded.regions$End))){
          msg <- 'options$excluded.regions should have integer End'
          stop(msg)
        }
      }else{
        msg <- 'Invalid options$excluded.regions. Please refer to ?ARTP2::options'
        stop(msg)
      }
    }
  }
  
  
  if(options$snp.miss.rate > .1 && !options$turn.off.filters){
    msg <- paste0("options$snp.miss.rate = ", options$snp.miss.rate, " might be too large")
    warning(msg)
  }
  
  if(options$inspect.snp.n <= 0){
    stop("options$inspect.snp.n should be a positive integer")
  }
  
  if(options$inspect.snp.percent < 0 || options$inspect.snp.percent > 1){
    stop("option$inspect.snp.percent should be in [0, 1]")
  }
  
  if(options$inspect.gene.n <= 0){
    stop("options$inspect.gene.n should be a positive integer")
  }
  
  if(options$inspect.gene.percent < 0 || options$inspect.gene.percent > 1){
    stop("option$inspect.gene.percent should be in [0, 1]")
  }
  
  if(length(intersect(options$excluded.snps, options$selected.snps)) > 0){
    stop("Some SNPs are specified in both options$excluded.snps and options$selected.snps")
  }
  
  if(!is.null(options$excluded.subs) && !is.null(options$selected.subs)){
    if(intersect(options$excluded.subs, options$selected.subs) > 0){
      stop("Some subject IDs are specified in both options$excluded.subs and options$selected.subs")
    }
  }
  
  if(options$gene.R2 <= 0){
    stop("gene.R2 should be in (0, 1]")
  }
  
  if(options$chr.R2 <= 0){
    stop("chr.R2 should be in (0, 1]")
  }
  
  if(options$huge.gene.R2 <= 0){
    stop("huge.gene.R2 should be in (0, 1]")
  }
  
  if(options$huge.chr.R2 <= 0){
    stop("huge.chr.R2 should be in (0, 1]")
  }
  
  
}
