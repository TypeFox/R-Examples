## COMMUNAL package
## part 1: run COMMUNAL across a range of data

#############  Main function, apply COMMUNAL over a range of data  #############

## v 1.1 added parallel functionality, and warnings if can't sort by variance
clusterRange <- function(dataMtx, ks, varRange, validation="all", verbose=T, ..., 
                         parallel=F, mc.cores=NULL, row.order=NULL) {
#   The default settings for COMMUNAL:
#   ks = 2:10
#   clus.methods = c("hierarchical", "kmeans"),
#   validation = c("Connectivity", "dunn", "wb.ratio", "g3", "g2", "pearsongamma", 
#                    "avg.silwidth", "sindex"),
#   dist.metric = "euclidean",
#   aggl.method = "average",
#   neighb.size = 10
  
  if (is.null(rownames(dataMtx))) {
    # add row names for indexing if necessary
    rownames(dataMtx) <- paste("Dim", 1:nrow(dataMtx), sep="")
    cat(paste("Added row names to dataMtx, from Dim1 to Dim", 
              nrow(dataMtx), "\n", sep=""))
  }
   
  if (is.null(row.order)) {
    # by default, sort by variance (keep highest-variance rows)
    data.var <- apply(dataMtx, 1, var, na.rm=TRUE)
  
    if(all(duplicated(data.var)[-1])){
      stop("Cannot sort by variances (all equal); please pass in ordering")
    } else if(any(duplicated(data.var))){
      warning("Some duplicated variances in object, sorting may not be optimal\nVariances:\t", 
              paste(data.var, collapse=", "), immediate.=TRUE)
    } 
    
    row.order <- order(data.var, decreasing=T)
  } 
  

  genes.list <- lapply(1:length(varRange), function(i) dataMtx[row.order[1:varRange[i]], , drop=FALSE])
  names(genes.list) <- paste0("vars_", unlist(varRange))
 
  ##  Run Function
  cat("Running COMMUNAL over range of variables...\t")
  if(parallel){
    if(is.null(mc.cores)) mc.cores <- parallel::detectCores()-1
    all.results <- parallel::mclapply(genes.list, mc.cores=mc.cores, mc.preschedule = FALSE,
                            mc.cleanup = FALSE, function(data) {
       cat("\n###################\n", nrow(data) ," variables \n")
       out <- COMMUNAL(data=data, parallel=T, mc.cores=mc.cores, ks=ks,
                       validation=validation, verbose=verbose, ...)
       return(out)
    })
  } else {
    all.results <- lapply(genes.list, function(data) {
      out <- COMMUNAL(data=data, parallel=F, ks=ks,
                      validation=validation, verbose=verbose, ...)
      cat("\n###################\n", nrow(data) ," variables complete \n")
      return(out)
    })
  }
  
  return(list(all.results=all.results, varRange=varRange))
}

