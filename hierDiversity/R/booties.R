booties <- 
function(dat, group, cluster, replace, reps, 
  q = 1, quant = c(0.025, 0.975), sims = FALSE) {
  if (is.matrix(group)) {
    Obs.div <- comp.div(dat, group, hier = (dim(group)[2] - 
      1), q = q, sims = FALSE)[[dim(group)[2]]]
    groupD <- data.frame(individual = 1:dim(group)[1], 
      group, dat)
    cluster <- c(cluster, "individual")
    replace <- c(replace, TRUE)
    temp.out <- matrix(NA, nrow = reps + 1, ncol = dim(Obs.div)[2])
    colnames(temp.out) <- names(Obs.div)
    temp.out <- data.frame(temp.out)
    temp.out[1, ] <- Obs.div
    temp.out[2:(reps + 1), ] <- do.call(rbind, 
      replicate(reps, expr = {
        temp <- reSample(dat = groupD, cluster = cluster, 
          replace = replace)
        DI <- dim(group)[2] + 1
        ngroup <- as.matrix(temp[, 2:DI])
        XX.div <- comp.div(dat = as.matrix(temp[, 
          -(1:DI)]), group = ngroup, hier = (DI - 
          2), sims = FALSE, q = q)
        XX.div[[length(cluster) - 1]][1, ]
      }, simplify = FALSE))
  }
  else {
    if (!is.null(dim(group))) {
      stop("something is wrong with grouping object")
    }
    Obs.div <- div.part(dat, group, q = q)
    temp.out <- matrix(NA, nrow = reps + 1, ncol = dim(Obs.div)[2])
    colnames(temp.out) <- names(Obs.div)
    temp.out <- data.frame(temp.out)
    temp.out[1, ] <- Obs.div
    for (i in 1:reps) {
      temp <- dat[sample(1:(dim(dat)[1]), dim(dat)[1], 
        replace = TRUE), ]
      temp.out[i + 1, ] <- div.part(dat = temp, 
        group = group, q = q)
    }
  }
  est <- matrix(NA, nrow = (length(quant) + 2), ncol = ncol(Obs.div[, 
    -1]))
  est <- data.frame(est, row.names = c(as.character(Obs.div[1, 
    1]), "SE", paste("q", quant, sep = "")))
  names(est) <- names((Obs.div[, -1]))
  est[1, ] <- Obs.div[, -1]
  est[2, ] <- apply(temp.out[, -1], 2, sd)
  for (i in 1:length(quant)) {
    est[i + 2, ] <- apply(temp.out[, -1], 2, quantile, 
      probs = quant[i], na.rm = TRUE)
  }
  if (sims) {
    list(estimates = est, boot.reps = temp.out)
  }
  else {
    return(est)
  }
}