
find.snps.in.regions <- function(stat, options){
  
  if(is.null(options$excluded.regions)){
    return(NULL)
  }
  
  exc.reg <- options$excluded.regions
  nfiles <- length(stat)
  exc.snps <- NULL
  comment <- NULL
  for(i in 1:nfiles){
    
    if(any(is.na(stat[[i]]$Chr))){
      msg <- 'Column \'Chr\' is missing or has NA in summary.files but options$excluded.regions is specified'
      stop(msg)
    }
    
    if(any(is.na(stat[[i]]$Pos))){
      msg <- 'Column \'Pos\' is missing or has NA in summary.files but options$excluded.regions is specified'
      stop(msg)
    }
    
    for(j in 1:nrow(exc.reg)){
      id <- which(stat[[i]]$Chr == exc.reg$Chr[j])
      if(length(id) == 0){
        next
      }
      st <- stat[[i]][id, c('SNP', 'Pos')]
      id <- which(st$Pos >= exc.reg$Start[j] & st$Pos <= exc.reg$End[j])
      if(length(id) == 0){
        next
      }
      exc.snps <- c(exc.snps, st$SNP[id])
      com <- paste0('Chr_', exc.reg$Chr[j], '_Start_', exc.reg$Start[j], '_End_', exc.reg$End[j])
      comment <- c(comment, rep(com, length(id)))
    }
  }
  
  if(!is.null(exc.snps)){
    tmp <- data.frame(exc.snps, comment, stringsAsFactors = FALSE)
    dup <- duplicated(tmp)
    tmp <- tmp[!dup, , drop = FALSE]
    exc.snps <- tmp$exc.snps
    comment <- tmp$comment
  }
  
  list(exc.snps = exc.snps, comment = comment)
  
}

