#### Code for jointIda (Preetam)
#### extract.parent.sets (version Nov 14)

extract.parent.sets <- function(x.pos, amat.cpdag, isCPDAG = FALSE) {
  amat.cpdag[which(amat.cpdag != 0)] <- 1
  amat.undir <- amat.cpdag*t(amat.cpdag)
  amat.dir <- amat.cpdag - amat.undir
  pasets.dir <- lapply(x.pos,function(x) which(amat.dir[,x] != 0))

  ## get all important connected components of the undirected subgraph
  conn.comp.imp <- NULL
  x.temp <- x.pos
  while (length(x.temp) > 0) {
    ## TODO: graph.dfs() -> dfs() {also in ../NAMESPACE !) once we rely on igraph >= 1.0.0
    comp.temp <- graph.dfs(graph = graph.adjacency(amat.undir, mode = 'undirected'),
                           root = x.temp[1], unreachable = FALSE)$order
    comp.temp <- comp.temp[!is.na(comp.temp)]
    x.temp <- setdiff(x.temp,comp.temp)
    conn.comp.imp <- c(conn.comp.imp, list(comp.temp))
  }
  ## new igraph 1.0.0 has more than just node numbers.  need only those
  ## "Hack":
  conn.comp.imp <- lapply(conn.comp.imp, as.integer)

  ## Chordality test, if required
  chordal <- if (!isCPDAG) {
    vapply(conn.comp.imp, function(i)
	   is.chordal(graph.adjacency(amat.undir[i,i], mode = "undirected"),
		      fillin = FALSE)$chordal,
	   NA)
  } else {
    rep(TRUE, length(conn.comp.imp))
  }

  ## Function for getting locally valid parent sets
  all.locally.valid.parents.undir <- function(amat,x) { # x must be a scaler
    amat.V <- as.integer(rownames(amat))
    pa.dir <- pasets.dir[[x.pos == amat.V[x]]]
    paset <- list(pa.dir)
    pa <- which(amat[,x] != 0) # cannot be a null set
    if (length(pa) == 1) {
      paset <- c(paset, list(c(amat.V[pa], pa.dir)))
    } else {
      for (i in 1:length(pa)) {
        pa.tmp <- combn(pa, i, simplify = TRUE)
        n.comb <- ncol(pa.tmp)
        for (j in 1:n.comb) {
          pa.t <- pa.tmp[, j]
          new.coll <-
            if (length(pa.t) > 1) {
              tmp <- amat[pa.t,pa.t]
              diag(tmp) <- 1
              (min(tmp) == 0)
            } else FALSE
          if (!new.coll)
            paset <- c(paset, list(c(amat.V[pa.t], pa.dir)))
        }
      }
    }
    return(paset)
  }

  extract.parent.sets.from.conn.comp <- function(i) {
    all.nodes <- conn.comp.imp[[i]]
    nvar <- length(all.nodes)
    if (nvar == 1) {
      pasets.comp <- list(pasets.dir[match(all.nodes,x.pos)])
    } else {
      conn.comp.mat <- amat.undir[all.nodes,all.nodes]
      rownames(conn.comp.mat) <- all.nodes
      x.. <- intersect(all.nodes, x.pos)
      m.x. <- match(x.., all.nodes)
      ii.x <- seq_along(x..) # = " 1:length(x..) "
      if(chordal[i] & nvar <= 12) {
        rownames(conn.comp.mat) <- colnames(conn.comp.mat) <- 1:nvar
        all.extensions <- allDags(conn.comp.mat, conn.comp.mat, NULL)
        pa.fun <- function(amat,j) c(all.nodes[which(amat[,m.x.[j]] != 0)],
                                     pasets.dir[[match(x..[j],x.pos)]])
        parent.sets.fun <- function(r) lapply(ii.x, pa.fun,
                                              amat = matrix(all.extensions[r,],nrow = nvar))
        pasets.comp <- lapply(1:nrow(all.extensions), parent.sets.fun)
      } else {
        pasets.comp <- lapply(m.x., all.locally.valid.parents.undir,
                              amat = conn.comp.mat)
        idx <- expand.grid(lapply(ii.x, function(j) seq_along(pasets.comp[[j]])))
        ## FIXME? Speed: interchange j and r
        pasets.comp <- lapply(1:nrow(idx), function(r)
                              lapply(ii.x, function(j) pasets.comp[[j]][[idx[r,j]]]))
      }

    }
    return(pasets.comp)
  }

  all.pasets <- lapply(1:length(conn.comp.imp), extract.parent.sets.from.conn.comp)
  idx <- expand.grid(lapply(1:length(all.pasets),
                            function(i) 1:length(all.pasets[[i]])))
  x.conn.comp <- unlist(lapply(1:length(conn.comp.imp),
                               function(i) intersect(conn.comp.imp[[i]], x.pos)))
  i.match <- match(x.pos, x.conn.comp)
  ## return all pasets :
  lapply(1:nrow(idx), function(i)
	 unlist(lapply(1:length(conn.comp.imp),
		       function(j) all.pasets[[j]][[idx[i,j]]]),
		recursive = FALSE)[i.match])
}## {extract.parent.sets}

## Main Function ----------------------------------------------------------
##
jointIda <- function(x.pos, y.pos, mcov, graphEst = NULL,
                     all.pasets = NULL, technique = c("RRC","MCD"))
{
  nx <- length(x.pos)
  if (is.null(all.pasets)) {
    amat <- as(graphEst,"matrix")
    amat[which(amat != 0)] <- 1
    all.pasets <- extract.parent.sets(x.pos,amat)
  } else { ## check format of all.pasets :
    if(!is.list(all.pasets) || vapply(all.pasets, length, 1L) != nx)
      stop("all.pasets is not given in an appropriate format.")
  }

  if (length(y.pos) > 1) { ## call myself for each y in y.pos :
    lapply(y.pos, function(y) jointIda(x.pos, y, mcov=mcov,
                                       all.pasets=all.pasets, technique=technique))
  } else {
    if (is.element(y.pos, x.pos))
      matrix(0, nrow = nx, ncol = length(all.pasets))
    else { ## return joint.effects
      technique <- match.arg(technique)
      switch(technique,
             RRC = matrix(unlist(lapply(all.pasets,function(pasets) RRC(mcov,x.pos,y.pos,pasets))),
                          nrow = nx),
             MCD = matrix(unlist(lapply(all.pasets,function(pasets) MCD(mcov,x.pos,y.pos,pasets))),
                          nrow = nx))
    }
  }
}

##' MCD :=  Modifying the Cholesky Decomposition
MCD <- function(cov.mat, intervention.set, var.y, pasets, return.modified.cov.mat = FALSE) {
  if (is.element(var.y,intervention.set) & !return.modified.cov.mat)
    return(rep(0,length(intervention.set)))
  if (length(intervention.set) == 1 & is.element(var.y,unlist(pasets)) & !return.modified.cov.mat)
    return(0)

  if (!return.modified.cov.mat) {
    imp.var <- unique(c(intervention.set,unlist(pasets),var.y))
    if (length(imp.var) < nrow(cov.mat)) {
      cov.mat <- cov.mat[imp.var,imp.var]
      intervention.set <- match(intervention.set,imp.var)
      var.y <- match(var.y,imp.var)
      pasets <- if(length(intervention.set) > 1)
        lapply(pasets, function(x) match(unlist(x), imp.var))
      else match(unlist(pasets), imp.var)
    }
  }

  do.Cholesky.modification <- function(x) {
    pa.x <- if(length(intervention.set) > 1) pasets[[match(x,intervention.set)]] else unlist(pasets)
    if (length(pa.x) == 0) return(cov.mat)
    ind <- c(pa.x,x,(1:nrow(cov.mat))[-c(pa.x,x)])
    cov.mat <- cov.mat[ind,ind]
    x <- match(x,ind)
    pa.x <- match(pa.x,ind)
    temp <- gchol(cov.mat)
    Lower.tri.mat <- solve(as.matrix(temp))
    tmp1 <- bdsmatrix::diag(temp)
    Diag.mat <- base::diag(tmp1)
    Lower.tri.mat[x,pa.x] <- 0
    ## MM{FIXME !} :
    cov.mat <- solve(Lower.tri.mat) %*% Diag.mat %*% t(solve(Lower.tri.mat))
    ## return
    cov.mat[order(ind),order(ind)]
  }

  for (i in 1:length(intervention.set)) {
    cov.mat <- do.Cholesky.modification(intervention.set[i])
  }

  if(return.modified.cov.mat)
    cov.mat
  else {
    MCD.estimate <- function(x) {
      if (is.element(var.y, unlist(pasets[match(x,intervention.set)])))
        0
      else
        cov.mat[var.y,x]/cov.mat[x,x]
    }
    vapply(intervention.set, MCD.estimate, double(1))
  }
} ## {MCD}


##' RRC := Recursive Regressions for Causal effects
RRC <- function(cov.mat, intervention.set, var.y, pasets) {

  adjusted.regression <- function(x,y) {
    if (x == y) return(0)
    pa.x <- if(length(intervention.set) > 1) pasets[[match(x,intervention.set)]] else unlist(pasets)
    if(is.element(y,pa.x)) 0 else
    solve(cov.mat[c(x,pa.x), c(x,pa.x)],
          cov.mat[c(x,pa.x), y])[1]
  }

  ## Define the vector of causal effects of intervention variables on var.y
  ## ISo   := intervention.set.on
  ## ISoIS := ISo.intervention.set := intervention.set.on.intervention.set
  ISo.var.y <- sapply(intervention.set, adjusted.regression, y = var.y)

  ## Define the required matrix of single intervention effects
  if (length(intervention.set) > 1) {
    ISoIS <-
      matrix(apply(expand.grid(intervention.set, intervention.set),
                   1L, function(x) adjusted.regression(x[1],x[2])),
             nrow = length(intervention.set))
  } else {
    return(ISo.var.y)
  }

  joint.effect.fun  <- function(x) {
    if(is.element(var.y, unlist(pasets[match(x, intervention.set)])))
      return(0)
    x.temp <- match(x, intervention.set)
    ## Intervention Set without x --> I.S. w/o x  --> IS.wo.x
    IS.wo.x <- intervention.set[-x.temp]
    ## Intitialize the RR estimate as the single intervention effect of intervention.set on var.y
    RR.estimate <- ISo.var.y[x.temp]
    ## Define the vector of causal effects of intervention.set on other intervention variables
    x.t.oIS <- ISoIS[x.temp,-x.temp]
    ## Define the vector of causal effects of "other" intervention variables on var.y
    ISo.var.y.temp <- ISo.var.y[-x.temp]
    ## Define the required matrix of single intervention effects
    ISoIS.temp <- ISoIS[-x.temp,-x.temp]

    while(length(IS.wo.x) > 1) {
      ## update RR.estimate and the other things accounting for
      ## the elimination of the first entry of the current intervention set
      RR.estimate <- RR.estimate - x.t.oIS[1]*ISo.var.y.temp[1]
      ISo.var.y.temp <- ISo.var.y.temp[-1] - ISoIS.temp[-1,1] * ISo.var.y.temp[1]
      x.t.oIS <- x.t.oIS[-1] - x.t.oIS[1]*ISoIS.temp[1,-1]
      if (length(IS.wo.x) > 2)
        ISoIS.temp <- ISoIS.temp[-1,-1] - tcrossprod(ISoIS.temp[-1,1],
                                                     ISoIS.temp[1,-1])
      IS.wo.x <- IS.wo.x[-1]
    }
    ## return
    RR.estimate - x.t.oIS * ISo.var.y.temp
  }

  sapply(intervention.set, joint.effect.fun)
}## {RRC}



###  MM: (ess-set-style 'DEFAULT) : we have much nesting ==> only indent by 2
## Local Variables:
## eval: (ess-set-style 'DEFAULT 'quiet)
## delete-old-versions: never
## End:

