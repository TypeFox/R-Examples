paint <- function (tree, subtree, branch, which = 1) {
  if (!is(tree,'ouchtree'))
    stop(sQuote("tree")," must be of class ",sQuote("ouchtree"))
  if (is(tree,'hansentree')) {
    regimes <- try(tree@regimes[[which]],silent=FALSE)
    if (inherits(regimes,'try-error'))
      stop(sQuote("paint")," error: invalid ",sQuote("which"))
  } else {
    regimes <- rep('unspec',length(tree@nodes))
    names(regimes) <- tree@nodes
  }
  if (!missing(subtree)) {
    st.nm <- names(subtree)
    if (is.null(st.nm))
      stop(sQuote("subtree")," must be a named vector")
    if (!all(st.nm%in%tree@nodes))
      stop("all names of ",sQuote("subtree")," must be names of nodes of ",sQuote("tree"))
    subtree <- as.character(subtree)
  } else {
    subtree <- character(0)
    st.nm <- character(0)
  }
  if (!missing(branch)) {
    br.nm <- names(branch)
    if (length(br.nm)>0) {
      if(is.null(br.nm))
        stop(sQuote("branch")," must be a named vector")
      if (!all(br.nm%in%tree@nodes))
        stop("all names of ",sQuote("branch")," must be names of nodes of ",sQuote("tree"))
      branch <- as.character(branch)
    }
  } else {
    branch <- character(0)
    br.nm <- character(0)
  }
  tog <- as.factor(c(as.character(subtree),as.character(branch)))
  subtree <- head(tog,length(subtree))
  branch <- tail(tog,length(branch))
  names(subtree) <- st.nm
  names(branch) <- br.nm
  for (k in seq(along=subtree)) {
    st <- sapply(tree@lineages,function(x,y)y%in%x[-1],y=st.nm[k])
    regimes[st] <- as.character(subtree[k])
  }
  for (k in seq(along=branch)) {
    br <- tree@nodes%in%br.nm[k]
    regimes[br] <- as.character(branch[k])
  }
  as.factor(regimes)
}

