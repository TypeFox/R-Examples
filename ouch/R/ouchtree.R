setClass(
         'ouchtree',
         representation=representation(
           nnodes = 'integer',
           nodes = 'character',
           ancestors = 'character',
           nodelabels = 'character',
           times = 'numeric',
           root = 'integer',
           nterm = 'integer',
           term = 'integer',
           anc.numbers = 'integer',
           lineages = 'list',
           epochs = 'list',
           branch.times = 'matrix',
           depth = 'numeric'
           )
         )

ouchtree <- function (nodes, ancestors, times, labels = as.character(nodes)) {

  nodes <- as.character(nodes)
  ancestors <- as.character(ancestors)

  n <- length(nodes)
  if (length(unique(nodes)) != n) stop("node names must be unique")
  if (length(ancestors) != n)
    stop("invalid tree: ",sQuote("ancestors")," must have the same length as ",sQuote("nodes"))
  if (length(times) != n) 
    stop("invalid tree: ",sQuote("times")," must have the same length as ",sQuote("nodes"))
  if (length(labels) != n)
    stop("invalid tree: ",sQuote("labels")," must be the same length as ",sQuote("nodes"))

  root <- which(is.root.node(ancestors))
  if (length(root) != 1)
    stop("invalid tree: there must be a unique root node, designated by its having ancestor = NA")
  if (times[root] != 0)
    stop("the algorithms assume that the root node is at time=0")
  
  term <- terminal.twigs(nodes,ancestors)
  if (length(term) <= 0) 
    stop("invalid tree: there ought to be at least one terminal node, don't you think?")

  outs <- which((!is.root.node(ancestors) & !(ancestors %in% nodes)))
  if (length(outs) > 0) {
    for (out in outs) {
      warning(
              sprintf(
                      "the ancestor of node %s is not in the tree",
                      nodes[out]
                      ),
              call.=FALSE)
    }
    stop("invalid tree")
  }
  
  anc <- ancestor.numbers(nodes,ancestors)

  if (any(anc==seq(along=anc),na.rm=TRUE)) {
    w <- which(anc==seq(along=anc))
    stop("this is no tree: node ",nodes[w[1]]," is its own ancestor",call.=FALSE)
  }

###  lineages <- build.lineages(anc)
  lineages <- vector(mode='list',length=n)
  todo <- root
  k <- 1
  while (k <= length(todo)) {
    todo <- c(todo,which(anc==todo[k]))
    a <- anc[todo[k]]
    lineages[[todo[k]]] <- c(todo[k],lineages[[a]])
    if (todo[k] %in% lineages[[a]]) 
      stop("this is no tree: circularity detected at node ",nodes[todo[k]]," in ",sQuote("ouchtree"),call.=FALSE)
    k <- k+1
  }

  for (k in 1:n) {
    if (!(root %in% lineages[[k]]))
      stop("node ",nodes[k]," is disconnected",call.=FALSE)
  }

  new(
      'ouchtree',
      nnodes = length(nodes),
      nodes=nodes,
      ancestors=ancestors,
      nodelabels=as.character(labels),
      times=times,
      root=root,
      nterm=length(term),
      term=term,
      anc.numbers=anc,
      lineages=lineages,
      epochs=epochs(lineages,times,term), # absolute times of epochs
      branch.times=branch.times(lineages,times,term), # absolute times of branch points
      depth=max(times)
      )
}

## map ancestor names to row numbers
ancestor.numbers <- function (nodenames, ancestors) { 
  sapply(ancestors,function(x)charmatch(x,nodenames),USE.NAMES=FALSE)
}

## nodenames of terminal twigs (terminal nodes are not ancestors)
terminal.twigs <- function (nodenames, ancestors) {
  which(nodenames %in% setdiff(nodenames,unique(ancestors)))
}

build.lineages <- function (ancestors) {
  n <- length(ancestors)
  lineages <- vector(mode='list',length=n)
  pedigree <- function (k) {
    if (is.null(lineages[[k]])) {
      a <- ancestors[k]
      if (is.root.node(a)) {
        lineages[[k]] <<- k
      } else {
        if (is.null(lineages[[a]])) Recall(a)
        if (k %in% lineages[[a]]) 
          stop('this is no tree: circularity detected at node ',k,call.=FALSE)
        lineages[[k]] <<- c(k,lineages[[a]])
      }
    }
    NULL
  }
  for (k in 1:n)
    pedigree(k)
  lineages
}

branch.times <- function (lineages, times, term) {
  N <- length(term)
  bt <- matrix(data=0,nrow=N,ncol=N)
  for (i in 1:N) {
    pedi <- lineages[[term[i]]]
    for (j in seq_len(i-1)) {
      pedj <- lineages[[term[j]]]
      for (k in 1:length(pedi)) {
        if (any(pedj == pedi[k])) break
      }
      bt[j,i] <- bt[i,j] <- times[pedi[k]]
    }
    bt[i,i] <- times[term[i]]
  }
  bt
}

epochs <- function (lineages, times, term) {
  N <- length(term)
  e <- vector(mode='list',length=N)
  for (i in 1:N) {
    pedi <- lineages[[term[i]]]
    e[[i]] <- times[pedi]
  }
  e
}

is.root.node <- function (anc) {
  is.na(anc)
}

