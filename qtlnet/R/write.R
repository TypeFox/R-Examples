###################################################################
## This creates object of class igraph.
write.qtlnet <- function(x,
                         filename,
                         edges = get.averaged.net(x, ...),
                         loci.list = loci.qtlnet(x, ...),
                         include.qtl=TRUE,
                         est.list = est.qtlnet(x, ...),
                         include.est=TRUE,
                         digits = 3,
                         col.names = TRUE,
                         ...)
{
  
  out <- cbind(cause = as.character(edges[[1]]),
               type = rep("Causal", nrow(edges)),
               effect = as.character(edges[[2]]))
  out.width <- edges[[3]]
  
  if(!is.null(loci.list) & include.qtl) {
    qtl <- unlist(loci.list)
    loci <- cbind(cause = qtl,
                  type = rep("QTL", length(qtl)),
                  effect = rep(names(loci.list), sapply(loci.list, length)))
    loci.width <- rep(1, nrow(loci))
    out <- rbind(out, loci)
    out.width <- c(out.width, loci.width)
  }
    
  out <- as.data.frame(out)
  out$width <- round(out.width, digits)

  if(include.est) {
    ## Need to figure out how to add estimates of edge coef.
    tmp <- rep(NA, nrow(out))
    for(effect in levels(out$effect)) {
      ii <- out$effect == effect
      if(any(ii)) {
        m <- match(as.character(out$cause)[ii], names(est.list[[effect]]),
                   nomatch = 0)
        tmp[ii][m > 0] <- est.list[[effect]][m]
      }
    }
    out$coef <- tmp
  }
  
  write.table(out, file = filename, quote = FALSE,
              col.names = col.names, row.names = FALSE)
  invisible(out)
}
