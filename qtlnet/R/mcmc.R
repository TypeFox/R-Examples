######################################################################
mcmc.qtlnet <- function(cross, pheno.col, threshold,
                        addcov=NULL, intcov=NULL,
                        nSamples = 1000, thinning=1,
                        max.parents = 3,
                        M0 = NULL,
                        burnin = 0.1, method = "hk", random.seed = NULL,
                        init.edges = 0,
                        saved.scores = NULL,
                        rev.method = c("nbhd", "node.edge", "single"),
                        verbose = FALSE, ...)
{
  ## Verbose: 1 or TRUE: saved count; 2: MCMC moves; 3: plot BIC; 4: 2&3.
  
  ## Check input parameters.
  rev.method <- match.arg(rev.method)

  ## Random number generator seed.
  if(!is.null(random.seed)) {
    if(!is.numeric(random.seed))
      stop("random seed must be numeric")
    set.seed(random.seed)
  }

  ## Burnin must be between 0 and 1.
  if(is.logical(burnin))
    burnin <- ifelse(burnin, 0.1, 0)
  if(burnin < 0 | burnin > 1)
    stop("burnin must be between 0 and 1")

  ## Adjust phenotypes and covariates to be numeric.
  cross <- adjust.pheno(cross, pheno.col, addcov, intcov)
  pheno.col <- cross$pheno.col
  pheno.names <- cross$pheno.names
  addcov <- cross$addcov
  intcov <- cross$intcov
  cross <- cross$cross

  n.pheno <- length(pheno.col)

  ## Initial network matrix.
  if(is.null(M0))
    M0 <- init.qtlnet(pheno.col, max.parents, init.edges)
  if(nrow(M0) != n.pheno | ncol(M0) != n.pheno)
    paste("M0 must be square matrix the size of pheno.col")

  ## LOD threshold by phenotype.
  if(length(threshold) == 1)
    threshold <- rep(threshold, n.pheno)
  if(length(threshold) != n.pheno)
    stop("threshold must have same length as pheno.col")

  ## Calculate genotype probabilities if not already done.
  if (!("prob" %in% names(cross$geno[[1]]))) {
    warning("First running calc.genoprob.")
    cross <- calc.genoprob(cross)
  }

  Mav <- matrix(0, n.pheno, n.pheno)
  n.burnin <- burnin * nSamples * thinning
  
  post.bic <- rep(NA,nSamples)
  post.model <- rep(NA,nSamples)
  all.bic <- rep(NA,nSamples)

  if(is.null(saved.scores))
    rev.method <- "single"
  saved.scores <- make.saved.scores(pheno.names, max.parents,
                                    saved.scores = saved.scores,
                                    verbose = verbose, ...)
        
  M.old <- M0
  ne.old <- nbhd.size(M.old, max.parents)[[1]]
  aux.new <- score.model(M.old, saved.scores, cross, addcov, intcov,
                         threshold, verbose, method = method, ...)
  
  bic.old <- aux.new$model.score
  tmp <- aux.new$update.scores
  codes <- dimnames(saved.scores)[[1]]
  if(!is.null(tmp)) {
    index <- nrow(saved.scores)
    index <- match(as.character(tmp$code), codes) + index * (tmp$pheno.col - 1)
    saved.scores[index] <- tmp$bic
  }
  model.old <- aux.new$model.name

  if(verbose) {
    cat("\n")
    if(verbose > 2)
      plot(c(1,nSamples), bic.old * c(0.5,1.2), type = "n",
           xlab = "sample", ylab = "BIC")
  }
  
  k <- 0
  if(thinning <= 1) {
    post.bic[1] <- all.bic[1] <- bic.old
    post.model[1] <- model.old
    k <- 1
    if(verbose > 2)
      points(k, post.bic[k], cex = 0.5)
  }
  cont.accept <- numeric(0)
  accept.fn <- function(x, move, fate = "accept") {
    move <- paste(fate, move, sep = ".")
    if(is.na(match(move, names(x))))
      x[move] <- 1
    else
      x[move] <- x[move] + 1
    x
  }
  for(i in 2:(nSamples*thinning)){
    ## Propose new network structure.
    M.new <- if(rev.method == "node.edge")
      propose.new.node.edge(M.old, max.parents,
                            saved.scores = saved.scores, rev.method = rev.method,
                            verbose = (verbose %in% c(2,4)), ...)
    else
      propose.new.structure(M.old, max.parents,
                            saved.scores = saved.scores, rev.method = rev.method,
                            verbose = (verbose %in% c(2,4)), ...)

    rev.ratio <- M.new$rev.ratio
    move <- M.new$move ## Want to keep track of moves and acceptance rates.
    M.new <- M.new$M
    
    ne.new <- nbhd.size(M.new, max.parents)[[1]]
    aux.new <- score.model(M.new, saved.scores, cross, addcov, intcov,
                           threshold, verbose, method = method, ...)

    ## Typically only a few update.scores.
    tmp <- aux.new$update.scores
    if(!is.null(tmp)) {
      index <- nrow(saved.scores)
      index <- match(as.character(tmp$code), codes) + index * (tmp$pheno.col - 1)
      saved.scores[index] <- tmp$bic
    }

    bic.new <- aux.new$model.score
    model.new <- aux.new$model.name

    ## Accept new model?
    ## Always if bic.new < bic.old.
    ## Otherwise do Metropolis-Hastings trick.
    if(is.infinite(rev.ratio))
      mr <- 1
    else
      mr <- exp(-0.5*(bic.new - bic.old)) * (ne.old / ne.new) * rev.ratio
    if(is.na(mr) | is.null(mr))
      browser()
    if(runif(1) <= min(1,mr)){
      M.old <- M.new
      bic.old <- bic.new
      ne.old <- ne.new
      model.old <- model.new
      cont.accept <- accept.fn(cont.accept, move)
    }
    else
      cont.accept <- accept.fn(cont.accept, move, "reject")

    ## Accumulate M's for average.
    if(i > n.burnin)
      Mav <- Mav + M.old
    
    ## Bookkeeping to save sample.
    if((i %% thinning) == 0){      
      k <- k + 1
      all.bic[k] <- bic.new
      post.bic[k] <- bic.old
      post.model[k] <- model.old
      if(verbose) {
        print(c(i,k)) 
        if(verbose > 2)
          points(k, post.bic[k], cex = 0.5)
      }
    }
  }

  ## Make average here.
  Mav <- Mav / (nSamples * thinning - n.burnin)

  cont.accept <- cont.accept[sort(names(cont.accept))]
  
  out <- list(post.model = post.model,
              post.bic = post.bic, 
              Mav = Mav,
              freq.accept = sum(cont.accept) / (nSamples * thinning),
              cont.accept = cont.accept,
              saved.scores=saved.scores,
              all.bic=all.bic,
              cross = cross)

  ## Attributes of qtlnet object.
  attr(out, "M0") <- M0
  attr(out, "threshold") <- threshold
  attr(out, "nSamples") <- nSamples
  attr(out, "thinning") <- thinning
  attr(out, "pheno.col") <- pheno.col
  attr(out, "pheno.names") <- pheno.names
  attr(out, "addcov") <- addcov
  attr(out, "intcov") <- intcov
  attr(out, "burnin") <- burnin
  attr(out, "method") <- method
  attr(out, "random.seed") <- random.seed
  attr(out, "random.kind") <- RNGkind()

  class(out) <- c("qtlnet","list")
  
  out
}
######################################################################
c.qtlnet <- function(...)
{
  ## Combine qtlnet objects.
  
  netlist <- list(...)
  out <- netlist[[1]]
  if(!inherits(out, "qtlnet")) {
    netlist <- out
    out <- netlist[[1]]
  }
  if(!inherits(out, "qtlnet"))
    stop("argument must be list of qtlnet objects")
  
  burnin <- attr(out, "burnin")
  n.pheno <- length(attr(out, "pheno.col"))
  
  if(length(netlist) > 1) {
    ## Need to do some checking here that qtlnet objects match up.
    ## Minimal for now.

    ## Assumes that parameters stay the same, although nSamples could change.

    if(n.pheno != mean(sapply(netlist, function(x) length(attr(x, "pheno.col")))))
      stop("different numbers of phenotypes not allowed")
    if(attr(out, "thinning") != mean(sapply(netlist, function(x) attr(x, "thinning"))))
      stop("thinning values differ")
    ## check threshold vector.
    if(!all(attr(out, "threshold") == apply(sapply(netlist,
                  function(x) attr(x, "threshold")), 1, mean)))
      stop("threshold values differ")
    if(!all(burnin == sapply(netlist, function(x) attr(x, "burnin"))))
      stop("burnin values differ")

    if(is.matrix(out$Mav)) {
      tmp <- out$Mav
      nsheet <- 1
      out$Mav <- array(0, c(nrow(tmp), ncol(tmp), length(netlist)))
    }
    else {
      ## Already is array. Want to add another sheet.
      tmp <- out$Mav
      dtmp <- dim(tmp)
      nsheet <- dtmp[3]
      dtmp[3] <- nsheet + length(netlist) - 1
      out$Mav <- array(0, dtmp)
    }
    out$Mav[, , seq(nsheet)] <- tmp
    
    out$cont.accept <- list(out$cont.accept)

    for(i in seq(2, length(netlist))) {
      ## Per step summaries.
      for(j in c("post.model","post.bic","all.bic","freq.accept"))
        out[[j]] <- c(out[[j]], netlist[[i]][[j]])

      out$cont.accept[[i]] <- netlist[[i]]$cont.accept
      
      ## Attributes.
      
      attr(out, "nSamples") <- c(attr(out, "nSamples"),
                                 attr(netlist[[i]], "nSamples"))

      ## Matrices of average network structures.
      out$Mav[ , , nsheet - 1 + i] <- netlist[[i]]$Mav
    }
  }
  out
}
######################################################################
nbhd.size <- function(M, max.parents = 3)
{
  n.deletions <- sum(M)
  n.additions <- 0 
  n.reversions <- 0 
  le <- ncol(M)
  for(j in 1:le){
        add.forbid <- forbidden.additions(M, j, max.parents)
        n.additions <- n.additions + (le - 1 - length(c(add.forbid$upf, add.forbid$downf)))
        rev.allow <- check.reversions(M, j, max.parents)
        if(!is.null(rev.allow))
          n.reversions <- n.reversions + nrow( rev.allow )
  }
  nbhd.size <- n.deletions + n.additions + n.reversions
  list(nbhd.size=nbhd.size, 
       n.deletions=n.deletions,
       n.additions=n.additions, 
       n.reversions=n.reversions)
}######################################################################
init.qtlnet <- function(pheno.col, max.parents = 3, init.edges = NULL)
{
  n.pheno <- length(pheno.col)
  if(is.null(init.edges)) {
    n.edges <- n.pheno * (n.pheno - 1) / 2
    init.edges <- sample(seq(0, n.edges), 1)
  }
  M <- matrix(0, n.pheno, n.pheno)
  if(init.edges > 0)
    for(i in seq(init.edges)) {
      cause <- sample(n.pheno, 1)
      if(i == 1)
        effect <- sample(seq(n.pheno)[-cause], 1)
      else
        effect <- propose.add(M, cause, max.parents)
      if(effect == 0)
        break
      M[cause, effect] <- 1
    }
  M
}
######################################################################
propose.add <- function(M, node, max.parents)
{
  ## Add an edge with a direction.
  aux1 <- forbidden.additions(M, node, max.parents)
  ## Is there a down option?
  aux3 <- unique(c(node, aux1$upf, aux1$downf))
  aux3 <- seq(ncol(M))[-aux3]
  
  if(length(aux3)) {
    ## Check if any of aux3 is at or above max.parents.

    wh <- which(apply(M[, aux3, drop = FALSE], 2, sum) >= max.parents)
    if(length(wh))
      aux3 <- aux3[-wh]
  }
  if(length(aux3) == 0)
    aux3 <- 0
  else {
    if(length(aux3) > 1)
      aux3 <- sample(aux3, 1)
  }
  aux3
}  
######################################################################
propose.new.node.edge <- function(M, max.parents = 3,
                                  saved.scores, rev.method = "node.edge",
                                  propose = c(node = 1, edge = 2, reverse = 10, drop = 2),
                                  check.down = FALSE,
                                  verbose = FALSE, ...)
{
  ## Ref: Grzegorczyk and Husmeier (2008) Mach Learn 71: 265–305.
  ## Extension is to updating individual nodes and dropping edges.
  rev.ratio <- 1

  ## Proposal weights for types of changes.
  if(any(propose <= 0))
    stop("propose values must be positive")
  name.propose <- names(propose)
  propose <- array(propose, 4)
  if(length(name.propose) == 4)
    names(propose) <- name.propose
  else
    names(propose) <- c("node","edge","reverse","drop")
  
  ## To speed up, use M as logical matrix. Still return 0/1 matrix for now.
  le.nodes <- ncol(M)
  le.edges <- sum(M)

  prob.node <- propose["node"] * le.nodes
  prob.node <- prob.node / (prob.node + propose["edge"] * le.edges)
  if(runif(1) <= prob.node) {
    move <- "node"
    ## Propose new parents for a node.
    node <- sample(seq(le.nodes), 1)
    if(verbose)
      cat(move, "")
    
    aux2 <- update.node(M, max.parents, saved.scores, node, verbose)
    rev.ratio <- aux2$rev.ratio
    M <- aux2$M
    if(verbose) {
      up <- node.parents(M, node)$parents
      cat(node, "up:", up, "\n")
    }
  }
  else {
    ## Propose change to edge.
    move <- "edge"
    if(verbose)
      cat(move, "")

    ## Pick a valid edge.
    edge <- which(M == 1)
    if(length(edge) > 1)
      edge <- sample(edge, 1)
    edge <- c(row(M)[edge], col(M)[edge])
    if(verbose)
      cat("edge ")
      
    ## Cyclic mistake checks.
    ## These should not happen!
    if(edge[1] == edge[2]) {
      cat("edge identity:", edge, "\n")
      browser()
    }
    aux1 <- check.downstream(M, edge[2])
    if(any(aux1 == edge[1])) {
      cat("edge downfall:", edge, aux1, "\n")
      browser()
    }
    
    prob.reverse <- propose["reverse"]
    prob.reverse <- prob.reverse / (prob.reverse + propose["drop"])
    if(runif(1) <= prob.reverse) {
      ## Propose to reverse an edge.
      if(verbose)  cat("reverse ")
      aux2 <- rev.edge(M, max.parents, saved.scores, edge, verbose)
    }
    else {
      ## Propose to drop an edge.
      if(verbose)  cat("drop ")
      aux2 <- drop.edge(M, max.parents, saved.scores, edge, verbose)
    }
    rev.ratio <- aux2$rev.ratio
    M <- aux2$M
    
    if(verbose) {
      for(i in 1:2) {
        up <- node.parents(M, edge[i])$parents
        cat(edge[i], "up:", up, "\n")
      }
    }

    if(check.down) {
      ## NOTE: This takes time and should be dropped eventually.
      ## Check if we somehow have a created a cycle earlier. This should not happen.
      down <- check.downstream(M, edge[2])[-1]
      up <- check.upstream(M, edge[2])[-1]
      if(length(down) & length(up)) {
        if(any(down %in% up)) {
          cat("downup:", edge, "\n", sort(down), "\n", sort(up), "\n")
          browser()
        }
      }
    }
  }

  
  list(M = M, rev.ratio = rev.ratio, move = move)
}
######################################################################
propose.new.structure <- function(M, max.parents = 3,
                                  saved.scores, rev.method = "nbhd",
                                  verbose = FALSE, ...)
{
  ## Acceptance rate is <20%. Could we improve here?
  ## Yes. See Grzegorczyk and Husmeier (2008) Mach Learn 71: 265–305.
  rev.ratio <- 1

  ## To speed up, use M as logical matrix. Still return 0/1 matrix for now.
  le.nodes <- ncol(M)
  flag <- TRUE
  while(flag){
    ## Pick node and decide on add/delete/reverse.
    ## Keep doing this until successful.
    
    node <- sample(seq(le.nodes),1)
    move <- sample(c("add","delete","reverse"),1)
    if(verbose)
      cat(node, move, "")

    switch(move,
           add = {
             aux3 <- propose.add(M, node, max.parents)
             if(!(flag <- (aux3 == 0))) {
               aux1 <- check.downstream(M, aux3)
               if(any(aux1 == node)) {
                 ## This should not happen!
                 cat(move, "downfall:", node, aux1, "\n")
                 browser()
               }
               M[node,aux3] <- 1
               if(verbose)
                 cat(aux3, "\n")
             } 
           },
           reverse = {
             ## Reverse direction of an edge.
             if(flag <- (rev.method == "single")) {
               aux1 <- check.reversions(M, node, max.parents)
               if(!(flag <- is.null(aux1))) {
                 le.rev <- nrow(aux1)
                 aux3 <- sample(seq(le.rev), 1)
                 aux3 <- aux1[aux3, ]
                 M[aux3[1], aux3[2]] <- 0
                 M[aux3[2], aux3[1]] <- 1
               }
             }
             else { ## Reverse method of Grzegorczyk and Husmeier.
               aux1 <- c(which(M[,node] == 1), le.nodes + which(M[node,] == 1))
               if(!(flag <- !length(aux1))) {
                 if(length(aux1) > 1)
                   aux1 <- sample(aux1, 1)
                 if(aux1 > le.nodes)
                   aux3 <- c(node, aux1 - le.nodes)
                 else
                   aux3 <- c(aux1, node)

                 aux1 <- check.downstream(M, aux3[2])
                 if(any(aux1 == aux3[1])) {
                   ## This should not happen!
                   cat(move, "downfall:", aux3[1], aux1, "\n")
                   browser()
                 }
                 
               }
               if(!flag) {
                 if(aux3[1] == aux3[2]) {
                   ## This should not happen!
                   cat(move, "identity:", aux3, "\n")
                   browser()
                 }
                 aux2 <- rev.edge(M, max.parents, saved.scores, aux3, verbose)
                 rev.ratio <- aux2$rev.ratio
                 M <- aux2$M
               }
             }
             if(!flag & verbose)
                 cat(aux3[1], aux3[2], "\n")
           },
           delete = {
             ## Delete an existing edge through node.
             aux1 <- c(which(M[, node] == 1), le.nodes + which(M[node,] == 1))
             if(!(flag <- !length(aux1))) {
               if(length(aux1) > 1)
                 aux1 <- sample(aux1, 1)
               if(aux1 > le.nodes)
                 aux3 <- c(node, aux1 - le.nodes)
               else
                 aux3 <- c(aux1, node)
               M[aux3[1], aux3[2]] <- 0
               if(verbose)
                 cat(aux3, "\n")
             }
           })
  }

  ## Check if we somehow have a created a cycle earlier. This should not happen.
  down <- check.downstream(M, aux3[2])[-1]
  up <- check.upstream(M, aux3[2])[-1]
  if(length(down) & length(up)) {
    if(any(down %in% up)) {
      cat(move, "downup:", aux3, "\n", sort(down), "\n", sort(up), "\n")
      browser()
    }
  }
  
  list(M = M, rev.ratio = rev.ratio, move = move)
}
######################################################################
update.node <- function(M, max.parents, saved.scores, node, verbose = FALSE)
{
  ## Scores for possible parent sets.
  down <- check.downstream(M, node)[-1]
  index <- index.parents(saved.scores, down, node)
  if(length(index))
    z.1 <- saved.scores[-index, node]
  else 
    z.1 <- saved.scores[, node]
  z.1 <- exp(min(z.1) - z.1)
  s.1 <- sum(z.1)
  if(length(z.1) == 0) ## Should not happen
    browser()
    
  ## Find probability for current parents of node.
  rev.ratio <- z.1[find.parent.score(M, saved.scores, -index, node)]
    
  ## Make node an orphan.
  M[,node] <- 0

  ## Sample new parents of node.
  parent <- sample(seq(length(z.1)), 1, prob = z.1)

  ## Proposal ratio.
  rev.ratio <- z.1[parent] / rev.ratio

  ## Add edges to graph for new parents of node.
  new.parent <- find.index.parent(M, saved.scores, -index, node, parent)
  if(length(new.parent))
    M[new.parent, node] <- 1

  if(verbose) {
    cat("node rev.ratio", rev.ratio, "\n")
    if(is.na(rev.ratio))
      browser()
  }
  
  list(M = M, rev.ratio = rev.ratio)

}
######################################################################
drop.edge <- function(M, max.parents, saved.scores, node.pair, verbose = FALSE)
{
  ## Grzegorczyk and Husmeier (2008) Mach Learn 71: 265–305.
  ## Call provides node.pair[1 -> 2].

  ## Step 1: Orphan both nodes after computing reverse proposal prob.
  rev.ratio <- reverse.proposal(M, saved.scores, node.pair, verbose)
  M[, node.pair] <- 0

  ## Step 2. Sample new parents for node 1 without node 2 as one parent.
  
  ## Exclude parents downstream of 1 or 2.
  tmp <- sample.exclude.down(M, saved.scores, node.pair, verbose)
  q.1 <- tmp$prob
  parent <- tmp$parent
  
  ## Add edges to graph for parents of node 1.
  if(length(parent))
    M[parent, node.pair[1]] <- 1

  ## Step 3. Sample new parents for node 2 without node 1 as parent.

  ## Exclude parents downstream of 1 or 2.
  tmp <- sample.exclude.down(M, saved.scores, rev(node.pair), verbose)
  q.2 <- tmp$prob
  parent <- tmp$parent

  ## Proposal ratio.
  rev.ratio <- q.1 * q.2 / rev.ratio

  ## Add edges to graph for parents of node 2.
  if(length(parent))
    M[parent, node.pair[2]] <- 1

  if(verbose) {
    cat("edge rev.ratio", rev.ratio, "\n")
    if(is.na(rev.ratio))
      browser()
  }
  
  list(M = M, rev.ratio = rev.ratio)
}
######################################################################
sample.exclude.down <- function(M, saved.scores, node.pair, verbose = FALSE)
{
  ## Exclude parents downstream of 1 or 2.
  down <- unique(c(check.downstream(M, node.pair[2]),
                   check.downstream(M, node.pair[1])[-1]))
  index <- index.parents(saved.scores, down, node.pair[1])
  z.1 <- saved.scores[-index, node.pair[1]]
  z.1 <- exp(min(z.1) - z.1)
  if(length(z.1) == 0) ## Should not happen
    browser()
  parent <- sample(seq(length(z.1)), 1, prob = z.1)
  q.1 <- z.1[parent] / sum(z.1)

  new.parent <- find.index.parent(M, saved.scores, -index, node.pair[1], parent)

  list(parent = new.parent, prob = q.1)
}
######################################################################
reverse.proposal <- function(M, saved.scores, node.pair, verbose = FALSE)
{
  ## Parents of 2 include 1 but exclude all downstream of 2.
  index <- index.parents(saved.scores, node.pair[1], node.pair[2])
  down <- check.downstream(M, node.pair[2])[-1]
  if(length(down))
    index <- index[-index.parents(saved.scores[index,, drop = FALSE], down, node.pair[2])]
  z.2 <- saved.scores[index, node.pair[2]]
  z.2 <- exp(min(z.2) - z.2)
  q.2 <- z.2[find.parent.score(M, saved.scores, index, node.pair[2])] / sum(z.2)

  ## Parents of 1 exclude nodes at or downstream of 2.
  down <- unique(c(check.downstream(M, node.pair[2]),
                   check.downstream(M, node.pair[1])[-1]))
  index <- index.parents(saved.scores, down, node.pair[1])
  z.1 <- saved.scores[-index, node.pair[1]]
  z.1 <- exp(min(z.1) - z.1)
  q.1 <- z.1[find.parent.score(M, saved.scores, -index, node.pair[1])] / sum(z.1)

  rev.ratio <- q.1 * q.2

  if(verbose) {
    cat("rev.ratio", rev.ratio, "\n")
    if(is.na(rev.ratio))
      browser()
  }
  rev.ratio
}
######################################################################
rev.edge <- function(M, max.parents, saved.scores, node.pair, verbose = FALSE)
{
  ## Grzegorczyk and Husmeier (2008) Mach Learn 71: 265–305.
  ## Call provides node.pair[1 -> 2].

  ## Step 1: Orphan both nodes after computing reverse proposal prob.
  rev.ratio <- reverse.proposal(M, saved.scores, node.pair, verbose)
  M[, node.pair] <- 0

  ## Step 2. Sample new parents for node 1 with node 2 as one parent.
  
  ## Exclude parents downstream of 1 (except 2).
  index <- index.parents(saved.scores, node.pair[2], node.pair[1])
  down <- check.downstream(M, node.pair[1])[-1]
  if(length(down))
    index <- index[-index.parents(saved.scores[index,, drop = FALSE], down, node.pair[1])]
  z.1 <- saved.scores[index, node.pair[1]]
  z.1 <- exp(min(z.1) - z.1)
  if(length(z.1) == 0) ## Should not happen
    browser()
  parent <- sample(seq(length(z.1)), 1, prob = z.1)
  q.1 <- z.1[parent] / sum(z.1)
  new.parent <- find.index.parent(M, saved.scores, index, node.pair[1], parent)

  ## Add edges to graph for parents of node 1.
  if(length(new.parent))
    M[new.parent, node.pair[1]] <- 1

  ## Step 3. Sample new parents for node 2 without node 1 as parent.

  ## Exclude parents downstream of 1 or 2.
  ## Somehow this is leading to cycles!
  ## SOmehow the rev(node.pair) is not working.
  tmp <- sample.exclude.down(M, saved.scores, rev(node.pair))
  q.2 <- tmp$prob
  new.parent <- tmp$parent

  ## Add edges to graph for parents of node 2.
  if(length(new.parent))
    M[new.parent, node.pair[2]] <- 1

  rev.ratio <- q.1 * q.2 / rev.ratio

  if(verbose) {
    cat("edge rev.ratio", rev.ratio, "\n")
    if(is.na(rev.ratio))
      browser()
  }
  
  list(M = M, rev.ratio = rev.ratio)
}
######################################################################
find.index.parent <- function(M, saved.scores, index, node, new.parent)
{
  if(length(index))
    parents <- dimnames(saved.scores)[[1]][index][new.parent]
  else
    parents <- dimnames(saved.scores)[[1]][new.parent]
  unlist(sapply(strsplit(parents, ",", fixed = TRUE),
                function(x, node) {
                  x <- as.numeric(x)
                  x + (x >= node)
                }, node))
}
######################################################################
index.parents <- function(saved.scores, node1, node2)
{
  ## Row names of saved.scores renumber around missing node.
  parents <- node1 - (node2 < node1)

  ## Which parents include node.pair[2]?
  which(sapply(strsplit(dimnames(saved.scores)[[1]], ",", fixed = TRUE),
               function(x, parents) any(x %in% parents),
               parents))
}
######################################################################
find.parent.score <- function(M, saved.scores, index, node)
{
  parents <- node.parents(M, node)$parents
  if(!is.null(parents)) {
    parents <- parents - (parents > node)
    if(length(index))
      match(paste(parents, collapse = ","), dimnames(saved.scores)[[1]][index])
    else
      match(paste(parents, collapse = ","), dimnames(saved.scores)[[1]])
  }
  else
    1
}
######################################################################
forbidden.additions <- function(M, node, max.parents = 3)
{
  ## upf are forbidden upstream nodes (already present downstream).
  ## downf are forbidden downstream nodes (already present upstream).
  ## Check on cycles is one step. See check.reversions() for longer cycles.

  ## Forbidden upstream additions
  upf <- which(M[node,] == 1)
  if(length(upf)==0)
    upf <- NULL

  ## Forbidden downstream additions. More complicated.
  downf <- which(M[,node] == 1)
  le <- length(downf)
  
  ## If at least max.parents are causal for node, forbid any more.
  if(le >= max.parents)
    downf <- seq(nrow(M))[-node]
  else {
    if(le > 0)
      downf <- check.upstream(M, downf)
    else
      downf <- NULL
  }
  
  list(upf=upf, downf=downf)
} 
######################################################################
check.qtlnet <- function(object,
                         min.prob = 0.9,
                         correct = TRUE,
                         verbose = FALSE,
                         ...)
{
  pheno.names <- attr(object, "pheno.names")
  n.pheno <- length(pheno.names)

  forbid <- 1
  while(!is.null(forbid)) {
    M <- threshold.net(object, min.prob = min.prob, ...)
    min.prob <- attr(M, "min.prob")
    M1 <- 1 * (M > 0)

    ## This may not be right yet.
    forbid <- NULL
    for(i in seq(n.pheno)) {
      downf <- check.upstream(M1, i)
      wh <- which(M1[i, downf] == 1)
      if(length(wh))
        forbid <- rbind(forbid, cbind(downf[wh],i, M[i, downf[wh]]))
    }
    if(!is.null(forbid)) {
      forbid <- data.frame(forbid)
      names(forbid) <- c("cause","react","prob")
      for(i in 1:2)
        forbid[[i]] <- ordered(pheno.names[forbid[[i]]], pheno.names)
    }
    if(verbose)
      print(forbid)
  }
  attr(M, "min.prob") <- min.prob
  
  list(forbid = forbid, M = M)
}
######################################################################
check.upstream <- function(M, nodes)
{
  ## After accounting for scanone, 85% of time is nbhd.size.
  ## Of that, 85% (75% overall) is in check.upstream.

  count.upok <- nrow(M) - length(nodes)
  
  ## check upstream to see if a cycle would be created.
  is.up <- apply(M[, nodes, drop = FALSE], 1, sum)
  flag <- TRUE
  while(flag){
    aux1 <- which(is.up > 0)
    if(flag <- (length(aux1) > 0 & count.upok > 0)) {
      aux1 <- aux1[!(aux1 %in% nodes)]
      if(flag <- length(aux1)) {
        is.up <- apply(M[, aux1, drop = FALSE], 1, sum)
        nodes <- c(nodes, aux1)
        count.upok <- count.upok - flag
      }
    }
  }
  nodes
}
######################################################################
check.downstream <- function(M, nodes)
{
  count.downok <- nrow(M) - length(nodes)
  
  ## check downstream to see if a cycle would be created.
  is.down <- apply(M[nodes,, drop = FALSE], 2, sum)
  flag <- TRUE
  while(flag){
    aux1 <- which(is.down > 0)
    if(flag <- (length(aux1) > 0 & count.downok > 0)) {
      aux1 <- aux1[!(aux1 %in% nodes)]
      if(flag <- length(aux1)) {
        is.down <- apply(M[aux1,, drop = FALSE], 2, sum)
        nodes <- c(nodes, aux1)
        count.downok <- count.downok - flag
      }
    }
  }
  nodes
}
######################################################################
check.reversions <- function(M, node, max.parents = 3)
{
  ## Check on possible reversals of directions.
  ## allowed if no cycles produced.
  ## forbidden if cycles result.

  allowed <- NULL
  
  ## Check upstream.
  up <- which(M[,node] == 1)
  le <- length(up)
  
  if(le) {
    if(le == 1)
      allowed <- cbind(up, node)
    else { ## le > 1
      ## Multiple colliders.
      forbid.up <- rep(FALSE, le)
      for(k in 1:le)
        forbid.up[k] <- up[k] %in% check.upstream(M, up[-k])
      if(any(forbid.up)){
        if(any(!forbid.up))
          allowed <- cbind(up[!forbid.up], node)
      }
      else
        allowed <- cbind(up, node)
    }
  }
    
  ## Check downstream.
  down <- which(M[node,] == 1)
  le <- length(down)
  
  allowed2 <- NULL
  
  if(le) {
    if(le == 1)
      allowed2 <- cbind(node, down)
    else { ## le > 1
      ## Multiple colliders.
      forbid.down <- rep(FALSE, le)
      for(k in 1:le)
        forbid.down[k] <- down[k] %in% check.downstream(M, down[-k])
      if(any(forbid.down)){
        if(any(!forbid.down))
          allowed2 <- cbind(node, down[!forbid.down])
      }
      else
        allowed2 <- cbind(down, node)
    }
  }
  allowed <- rbind(allowed, allowed2)

  ## Final check of max.parents.
  if(!is.null(allowed)) {
    wh <- which(apply(M[, allowed[,1], drop = FALSE], 2, sum) >= max.parents)
    if(length(wh)) {
      ## Remove pairs.
      allowed <- allowed[-wh,, drop = FALSE]
      if(nrow(allowed) == 0)
        allowed <- NULL
    }
  }
    
  allowed
}
##########################################################################
legacy.qtlnet <- function(qtlnet.object, codes = TRUE)
{
  ## Convert old qtlnet object to current format. Idempotent.
  
  qtlnet.object$Mav <- get.model.average(qtlnet.object)
  qtlnet.object$post.net.str <- NULL
  qtlnet.object$model.average <- NULL
  if(codes)
    dimnames(qtlnet.object$saved.scores)[[1]] <-
      legacy.code(dimnames(qtlnet.object$saved.scores)[[1]])
  qtlnet.object
}
