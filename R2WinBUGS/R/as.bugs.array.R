as.bugs.array <- function(sims.array, model.file=NULL, program=NULL,
                          DIC=FALSE, DICOutput=NULL,
                          n.iter=NULL, n.burnin=0, n.thin=1)
{
  ## Jouni Kerman's function to convert a 3-way array to a Bugs object
  ##
  ## 'sims.array' is supposed to be a 3-way array with
  ## n.sims*n.chains*n.parameters simulations, and
  ## the 3rd component of dimnames(x) should have the parameter names.
  d <- dim(sims.array)
  n.keep       <- d[1]
  n.chains     <- d[2]
  n.parameters <- d[3]
  n.sims       <- n.keep*n.chains
  if (is.null(n.iter)) n.iter <- (n.keep + n.burnin)*n.thin
  #
  parameter.names <- dimnames(sims.array)[[3]]
  if (is.null(parameter.names)) {
    parameter.names <- paste("P", 1:n.parameters, sep="")
    dimnames(sims.array)[[3]] <- parameter.names
  }
  parameters.to.save <- unique(sapply(strsplit(parameter.names,"\\["),"[",1))
  #
  sims <- matrix(NA, n.sims, n.parameters)
  root.long <- character(n.parameters)
  indexes.long <- vector(n.parameters, mode = "list")
  for (i in 1:n.parameters) {
    temp <- decode.parameter.name(parameter.names[i])
    root.long[i] <- temp$root
    indexes.long[[i]] <- temp$indexes
  }
  n.roots <- length(parameters.to.save)
  left.bracket.short <- as.vector(regexpr("[[]", parameters.to.save))
  right.bracket.short <- as.vector(regexpr("[]]", parameters.to.save))
  root.short <- ifelse(left.bracket.short == -1, parameters.to.save,
    substring(parameters.to.save, 1, left.bracket.short - 1))
  dimension.short <- rep(0, n.roots)
  indexes.short <- vector(n.roots, mode = "list")
  n.indexes.short <- vector(n.roots, mode = "list")
  long.short <- vector(n.roots, mode = "list")
  length.short <- numeric(n.roots)
  for (j in 1:n.roots) {
    long.short[[j]] <- (1:n.parameters)[root.long == root.short[j]]
    length.short[j] <- length(long.short[[j]])
    if (length.short[j] == 0){
      stop(paste("parameter", root.short[[j]], "is not in the model"))
    }
    else if (length.short[j] > 1) {
      dimension.short[j] <- length(indexes.long[[long.short[[j]][1]]])
      n.indexes.short[[j]] <- numeric(dimension.short[j])
      for (k in 1:dimension.short[j]){
        n.indexes.short[[j]][k] <- length(unique(unlist(lapply(indexes.long[long.short[[j]]], .subset, k))))
      }
      length.short[j] <- prod(n.indexes.short[[j]])
### bug reported by S. McKay Curtis on February 22, 2010, we cannot check that safely:
##      if (length(long.short[[j]]) != length.short[j]){
##        stop(paste("error in parameter", root.short[[j]],
##          "in parameters.to.save"))
##      }
      indexes.short[[j]] <- as.list(numeric(length.short[j]))
      for (k in 1:length.short[j]){
        indexes.short[[j]][[k]] <- indexes.long[[long.short[[j]][k]]]
      }
    }
  }
  rank.long <- unlist(long.short)
  # -----
  # yes, it's inefficient to do this, but for now I'm just letting this be as it is:
  for (k in 1:n.parameters) {
    sims[,k] <- as.vector(sims.array[,,k])
  }
  # ----
  dimnames(sims) <- list(NULL, parameter.names)
  summary <- monitor(sims.array, n.chains, keep.all = TRUE)
  last.values <- as.list(numeric(n.chains))
  for (i in 1:n.chains) {
    n.roots.0 <- if(DIC) n.roots - 1 else n.roots
    last.values[[i]] <- as.list(numeric(n.roots.0))
    names(last.values[[i]]) <- root.short[1:n.roots.0]
    for (j in 1:n.roots.0) {
      if (dimension.short[j] <= 1) {
        last.values[[i]][[j]] <- sims.array[n.keep, i, long.short[[j]]]
        names(last.values[[i]][[j]]) <- NULL
      }
      else{
        last.values[[i]][[j]] <- aperm(array(sims.array[n.keep, i,
          long.short[[j]]], rev(n.indexes.short[[j]])), dimension.short[j]:1)
      }
    }
  }
  sims <- sims[sample(n.sims), , drop = FALSE]
  sims.list <- summary.mean <- summary.sd <- summary.median <- vector(n.roots,
    mode = "list")
  names(sims.list) <- names(summary.mean) <- names(summary.sd) <- names(summary.median) <- root.short
  for (j in 1:n.roots) {
    if (length.short[j] == 1) {
      sims.list[[j]] <- sims[, long.short[[j]]]
      summary.mean[[j]] <- summary[long.short[[j]], "mean"]
      summary.sd[[j]] <- summary[long.short[[j]], "sd"]
      summary.median[[j]] <- summary[long.short[[j]], "50%"]
    }
    else {
      temp2 <- dimension.short[j]:1
      sims.list[[j]] <- aperm(array(sims[, long.short[[j]]],
        c(n.sims, rev(n.indexes.short[[j]]))), c(1, (dimension.short[j]+1):2))
      summary.mean[[j]] <- aperm(array(summary[long.short[[j]],
        "mean"], rev(n.indexes.short[[j]])), temp2)
      summary.sd[[j]] <- aperm(array(summary[long.short[[j]],
        "sd"], rev(n.indexes.short[[j]])), temp2)
      summary.median[[j]] <- aperm(array(summary[long.short[[j]],
        "50%"], rev(n.indexes.short[[j]])), temp2)
    }
  }
  summary <- summary[rank.long, ]

  all <- list(n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
    n.thin = n.thin, n.keep = n.keep, n.sims = n.sims,
    sims.array = sims.array[,,rank.long,drop = FALSE], sims.list = sims.list,
    sims.matrix = sims[, rank.long], summary = summary, mean = summary.mean,
    sd = summary.sd, median = summary.median, root.short = root.short,
    long.short = long.short, dimension.short = dimension.short,
    indexes.short = indexes.short, last.values = last.values, program=program,
    model.file=model.file)

  if(DIC && is.null(DICOutput)) { ## calculate DIC from deviance
    deviance <- all$sims.array[, , dim(sims.array)[3], drop = FALSE]
    dim(deviance) <- dim(deviance)[1:2]
    pD <- numeric(n.chains)
    DIC <- numeric(n.chains)
    for (i in 1:n.chains) {
      pD[i] <- var(deviance[, i])/2
      DIC[i] <- mean(deviance[, i]) + pD[i]
    }
    all <- c(all, list(isDIC=TRUE, DICbyR=TRUE,
                       pD=mean(pD), DIC=mean(DIC)))
  } else if(DIC && !is.null(DICOutput)) { ## use DIC from BUGS
    all <- c(all, list(isDIC=TRUE, DICbyR=FALSE,
                       pD=DICOutput[nrow(DICOutput),4],
                       DIC=DICOutput[nrow(DICOutput),3]))
  } else {
    all <- c(all, isDIC=FALSE)
  }
  class(all) <- "bugs"
  all
}
