flimSD <-
function(fbo, response, grouping=NULL) {
  samples <- fbo$samples
  bml <- lapply(samples, flimMean, response = response,
                grouping = grouping)
  bmm <- do.call("cbind", bml)
  if(is.null(grouping)) {
    bmmh <- bmm[, seq(1, dim(bmm)[2]-1, 2)]
    stds <- apply(bmmh, 1, sd)
    stds
  } else {
    g.levels <- unique((fbo$org$df)[, grouping])
    bmmh <- bmm[, names(bmm)%in%g.levels]
    no.samples <- length(samples)
    nr.col <- no.samples*length(g.levels)
    sd.list <- list()
    for(l in 1:length(g.levels)) {
      idx <- seq(l, nr.col, length(g.levels))
      #print(bmmh[, idx])
      bmmh.l.sd <- apply(bmmh[, idx], 1, sd)
      sd.list[[l]] <- bmmh.l.sd
    }
    stds <- as.data.frame(do.call("cbind", sd.list))
    names(stds) <- g.levels
    stds
  }
}
