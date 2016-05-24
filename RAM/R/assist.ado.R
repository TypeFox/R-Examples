assist.ado <- function(data, meta, is.OTU=TRUE, ranks=NULL,           
                       data.trans=NULL, dist=NULL, meta.strata=NULL, 
                       perm=1000, top=NULL, mode="number") {
  if ( is.OTU ) {
    .valid.data(data, is.OTU=is.OTU)
  }
  num.data <- length(data)
  labels <- names(data)
  matched <- match.data(data, is.OTU=is.OTU, meta=meta)
  
  # rework with the data based on OTU or tax.abund matrix
  data.new <- data.revamp(data=data, is.OTU=is.OTU, ranks=ranks, 
                          stand.method=data.trans, top=top, mode=mode)
  
  
  # filter METAdata, exclude variable with only 1 level, with
  # missing data, and non numeric&&factor/charactor (NNF)
  suppressWarnings(meta.new <- filter.META(meta, 
                                           exclude=meta.strata))
  
  if ( is.null(meta.strata) ) {
    strata <- NULL
  } else {
    strata <- meta.new[[meta.strata]]
  }
  
  
  # adonis
  # either use distance matrix of the ecology data
  # or use the ecology data with or without being standardized.
  ado.list <- list()
  for ( i in 1:length(data.new) ) {
    dt <- data.new[[i]]
    if ( is.null(dt) ) { break }
    label <- names(data.new)[i]
    
    if ( is.null(dist) ) {
      ado <- vegan::adonis(dt ~ ., data=meta.new, 
                           permutations=perm, 
                           strata= strata)
    } else {
      ado <- vegan::adonis(dt ~ ., data=meta.new, 
                           method=dist, 
                           permutations=perm, 
                           strata=strata) 
    }
    #names(data.ado)[i] <- names(data.new)[i]
    ado.list[[label]] <- ado
  }
  
  return(ado.list)
}
