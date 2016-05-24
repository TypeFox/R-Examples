cat.pdb <- function(..., renumber=FALSE, rechain=TRUE) {
  cl <- match.call()
   
  objs <- list(...)
  are.null <- unlist(lapply(objs, is.null))
  objs <- objs[!are.null]

  if(length(objs)==1) 
     if(is.null(cl$rechain)) rechain = FALSE
  else
     if(length(objs)<1) return(NULL)

  if(any(!unlist(lapply(objs, is.pdb))))
    stop("provide PDB objects as obtained from read.pdb()")
  
  ## avoid NA as chain ID
  na.inds <- lapply(objs, function(x) is.na(x$atom$chain))
  if(any(unlist(na.inds))) {
    na.inds <- which(unlist(lapply(na.inds, function(x) any(x))))
    for(i in na.inds) {
      tmp <- objs[[i]]
      tmp$atom$chain[ is.na(tmp$atom$chain) ] <- " "
      objs[[i]] <- tmp
    }
  }

  ## save original chain IDs
  ori.chain <- unlist(lapply(objs, function(x) unique(x$atom$chain)))

  ## always assign new chain identifiers 
  ## and bring back original chain ID later if rechain=FALSE
  k <- 1
  chain.repo <- c(LETTERS, letters, 0:9)
  for(i in 1:length(objs)) {
    x <- objs[[i]]
    objs[[i]] <- .update.chain(x, chain.repo[k:length(chain.repo)])
    k <- k + length(unique(x$atom$chain))
  }

  ## concat objects
  new <- objs[[1]]
  if(length(objs) > 1) { 
     for(i in 2:length(objs)) {
       new$atom     <- rbind(new$atom, objs[[i]]$atom)
       new$xyz      <- cbind(new$xyz, objs[[i]]$xyz)
       new$seqres <- c(new$seqres, objs[[i]]$seqres)
       new$helix <- c(new$helix, objs[[i]]$helix)
       new$sheet <- c(new$sheet, objs[[i]]$sheet)
     }
  }
  ## merge SSE info
  for(i in c("helix", "sheet")) {
     sse <- new[[i]]
     if(!is.null(sse)) {
        coms <- names(sse)
        names(sse) <- NULL # avoid nested naming in results
        inds <- which(!duplicated(coms))
        for(j in inds) 
           sse[[j]] <- do.call(c, sse[coms %in% coms[j]])
        sse <- sse[inds]
        names(sse) <- coms[inds]
        new[[i]] <- sse
     }
  }

  ## renumber residues
  chk <- try(clean.pdb(new, consecutive = !rechain, 
                force.renumber = renumber, verbose=FALSE), silent=TRUE)
  if(inherits(chk, "try-error")) {
     warning("cat.pdb(): Bad format pdb generated. Try rechain=TRUE and/or renumber=TRUE")
     new['helix'] <- list(NULL)
     new['sheet'] <- list(NULL)
  }
  else
     new <- chk

  if(!rechain) new <- .update.chain(new, ori.chain)

  ## build new PDB object
  new$call <- cl

  ## remap " " chain IDs to NA values
  new$atom$chain[ new$atom$chain==" " ] <- as.character(NA)

  ## check connectivity
  chains <- unique(new$atom$chain)
  for(i in 1:length(chains)) {
    sele <- atom.select(new, chain=chains[i], verbose=FALSE)
    tmp <- trim.pdb(new, sele)
    
    if(!inspect.connectivity(tmp))
      warning(paste("possible chain break in molecule: chain", chains[i]))
  }
  
  return(new)
}

.update.chain <- function(x, chain.repo = LETTERS) {
   new <- x
   chains <- unique(x$atom$chain)
   for(j in 1:length(chains)) {
     inds <- which(x$atom$chain==chains[j])
     new$atom$chain[inds] <- chain.repo[j]
     if(!is.null(x$helix)) new$helix$chain[x$helix$chain==chains[j]] <- chain.repo[j]
     if(!is.null(x$sheet)) new$sheet$chain[x$sheet$chain==chains[j]] <- chain.repo[j]
   }
   new
}
