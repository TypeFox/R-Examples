as.party <- function(obj, ...)
    UseMethod("as.party")

as.party.rpart <- function(obj, ...) {

    ff <- obj$frame
    n  <- nrow(ff)
    if (n==1) return(partynode(as.integer(1)))  # special case of no splits

    is.leaf <- (ff$var == "<leaf>")
    vnames <- ff$var[!is.leaf]  #the variable names for the primary splits

    index <- cumsum(c(1, ff$ncompete + ff$nsurrogate + 1*(!is.leaf)))
    splitindex <- list()
    splitindex$primary <- numeric(n)
    splitindex$primary[!is.leaf] <- index[c(!is.leaf, FALSE)]
    splitindex$surrogate <- lapply(1L:n, function(i) {
        prim <- splitindex$primary[i]
        if (prim < 1 || ff[i, "nsurrogate"] == 0) return(NULL)
        else return(prim + ff[i, "ncompete"] + 1L:ff[i, "nsurrogate"])
    })
    
    mf <- model.frame(obj)
    
    ## check if any of the variables in the model frame is a "character"
    ## and convert to "factor" if necessary
    for(i in which(sapply(mf, function(x) class(x)[1L]) == "character")) mf[[i]] <- factor(mf[[i]])

    rpart_fitted <- function() {
	ret <- as.data.frame(matrix(nrow = NROW(mf), ncol = 0))
	ret[["(fitted)"]] <- obj$where
        ret[["(response)"]] <- model.response(mf)        
        ret[["(weights)"]] <- model.weights(mf)
        ret
    }
    fitted <- rpart_fitted()

    rpart_kids <- function(i) {
        if (is.leaf[i]) return(NULL)
        else return(c(i + 1L, 
            which((cumsum(!is.leaf[-(1L:i)]) + 1L) == cumsum(is.leaf[-(1L:i)]))[1L] + 1L + i))
    }

    rpart_onesplit <- function(j) {
        if (j < 1) return(NULL)
        ### numeric
        if (abs(obj$split[j, "ncat"]) == 1) {
            ret <- partysplit(varid = which(rownames(obj$split)[j] == names(mf)),
                      breaks = as.double(obj$split[j, "index"]),
                      right = FALSE,
                      index = if(obj$split[j, "ncat"] > 0) 2:1)
        } else {
            index <- obj$csplit[obj$split[j, "index"],]
            ### csplit has columns 1L:max(nlevels) for all factors
            ### index <- index[1L:obj$split[j, "ncat"]] ??? safer ???
            index <- index[1L:nlevels(mf[, rownames(obj$split)[j]])]
            index[index == 2L] <- NA ### level not present in split
            index[index == 3L] <- 2L  ### 1..left, 3..right
            ret <- partysplit(varid = which(rownames(obj$split)[j] == names(mf)),
                      index = as.integer(index))
        }
        ret
    }
                      
    rpart_split <- function(i)
        rpart_onesplit(splitindex$primary[i])
    
    rpart_surrogates <- function(i)
        lapply(splitindex$surrogate[[i]], rpart_onesplit)

    rpart_node <- function(i) {
        if (is.null(rpart_kids(i))) return(partynode(as.integer(i)))
        nd <- partynode(as.integer(i), split = rpart_split(i),
	           kids = lapply(rpart_kids(i), rpart_node),
	           surrogates = rpart_surrogates(i))

        ### determine majority for (non-random) splitting
        left <- nodeids(kids_node(nd)[[1L]], terminal = TRUE)
        right <- nodeids(kids_node(nd)[[2L]], terminal = TRUE)
        nd$split$prob <- c(0, 0)
        nl <- sum(fitted[["(fitted)"]] %in% left)
        nr <- sum(fitted[["(fitted)"]] %in% right)
        nd$split$prob <- if (nl > nr) c(1, 0) else c(0, 1)
        nd$split$prob <- as.double(nd$split$prob)
        return(nd)
    }

    node <- rpart_node(1)

    rval <- party(node = node, data = mf[0L,], fitted = fitted,
      terms = obj$terms, info = list(method = "rpart"))
    class(rval) <- c("constparty", class(rval))
    return(rval)
}

model.frame.rpart <- function(formula, ...) {
  ## if model.frame is stored, simply extract
  if(!is.null(formula$model)) return(formula$model)
  
  ## otherwise reevaluate model.frame using original call
  mf <- formula$call
  mf <- mf[c(1L, match(c("formula", "data", "subset", "na.action", "weights"), names(mf), 0L))]
  if (is.null(mf$na.action)) mf$na.action <- rpart::na.rpart
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  
  ## use terms instead of formula in call
  mf$formula <- formula$terms
  
  ## evaluate in the right environment and return
  env <- if(!is.null(environment(formula$terms))) environment(formula$terms) else parent.frame()
  mf <- eval(mf, env)
  return(mf)
}

as.party.Weka_tree <- function(obj, ...) {

  ## needs RWeka and rJava
  stopifnot(requireNamespace("RWeka"))

  ## J48 tree? (can be transformed to "constparty")
  j48 <- inherits(obj, "J48")

  ## construct metadata
  mf <- model.frame(obj)
  mf_class <- sapply(mf, function(x) class(x)[1L])
  mf_levels <- lapply(mf, levels)

  x <- rJava::.jcall(obj$classifier, "S", "graph")
  if(j48) {
    info <- NULL
  } else {
    info <- RWeka::parse_Weka_digraph(x, plainleaf = FALSE)$nodes[, 2L]
    info <- strsplit(info, " (", fixed = TRUE)
    info <- lapply(info, function(x) if(length(x) == 1L) x else c(x[1L], paste("(", x[-1L], sep = "")))
  }
  x <- RWeka::parse_Weka_digraph(x, plainleaf = TRUE)
  nodes <- x$nodes
  edges <- x$edges
  is.leaf <- x$nodes[, "splitvar"] == ""

  weka_tree_kids <- function(i) {
    if (is.leaf[i]) return(NULL)
      else return(which(nodes[,"name"] %in% edges[nodes[i,"name"] == edges[,"from"], "to"]))
  }

  weka_tree_split <- function(i) {
    if(is.leaf[i]) return(NULL)
    
    var_id <- which(nodes[i, "splitvar"] == names(mf))
    edges <- edges[nodes[i,"name"] == edges[,"from"], "label"]
    split <- Map(c, sub("^([[:punct:]]+).*$", "\\1", edges), sub("^([[:punct:]]+) *", "", edges))
    ## ## for J48 the following suffices
    ## split <- strsplit(edges[nodes[i,"name"] == edges[,"from"], "label"], " ")

    if(mf_class[var_id] %in% c("ordered", "factor")) {
      stopifnot(all(sapply(split, head, 1) == "="))
      stopifnot(all(sapply(split, tail, 1) %in% mf_levels[[var_id]]))
      
      split <- partysplit(varid = as.integer(var_id),
        index = match(mf_levels[[var_id]], sapply(split, tail, 1)))
    } else {
      breaks <- unique(as.numeric(sapply(split, tail, 1)))
      breaks <- if(mf_class[var_id] == "integer") as.integer(breaks) else as.double(breaks) ## FIXME: check
      
      stopifnot(length(breaks) == 1 && !is.na(breaks))
      stopifnot(all(sapply(split, head, 1) %in% c("<=", ">")))
      
      split <- partysplit(varid = as.integer(var_id),
        breaks = breaks, right = TRUE,
	index = if(split[[1L]][1L] == ">") 2L:1L)
    }
    return(split)
  }

  weka_tree_node <- function(i) {
    if(is.null(weka_tree_kids(i))) return(partynode(as.integer(i), info = info[[i]]))
    partynode(as.integer(i),
      split = weka_tree_split(i),
      kids = lapply(weka_tree_kids(i), weka_tree_node))
  }

  node <- weka_tree_node(1)

  if(j48) {
    pty <- party(
      node = node,
      data = mf[0L,],
      fitted = data.frame("(fitted)" = fitted_node(node, mf),
        		  "(response)" = model.response(mf),
        		  check.names = FALSE),
      terms = obj$terms,
      info = list(method = "J4.8"))
    class(pty) <- c("constparty", class(pty))
  } else {
    pty <- party(
      node = node,
      data = mf[0L,],
      fitted = data.frame("(fitted)" = fitted_node(node, mf), check.names = FALSE),
      terms = obj$terms,
      info = list(method = class(obj)[1L]))      
  }

  return(pty)
}
