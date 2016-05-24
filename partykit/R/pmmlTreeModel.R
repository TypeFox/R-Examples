pmmlTreeModel <- function(file, ...) {
  stopifnot(requireNamespace("XML"))
  as.party(XML::xmlRoot(XML::xmlTreeParse(file)))
}

as.party.XMLNode <- function(obj, ...) {

  stopifnot(requireNamespace("XML"))
  ## check whether XML specifies a TreeModel
  stopifnot(c("DataDictionary", "TreeModel") %in% names(obj))
  if(any(warnx <- c("MiningBuildTask", "TransformationDictionary", "Extension") %in% names(obj)))
    warning(sprintf("%s not yet implemented", paste(names(obj)[warnx], collapse = ", ")))
  
  ## process header information
  if("Header" %in% names(obj)) {
    hdr <- obj[["Header"]]
    h_info <- c(Header = paste(as.character(XML::xmlAttrs(hdr)), collapse = ", "))
    if(length(hdr) > 0L) {
      h_info <- c(h_info,
        XML::xmlSApply(hdr, function(x)
	  paste(c(as.character(XML::xmlAttrs(x)), XML::xmlValue(x)), collapse = ", ")
	)
      )
    }
  } else {
    h_info <- NULL
  }
  
  ## parse data dictionary
  extract_empty_model_frame <- function(x) {

    ## extract DataDictionary
    dd <- x[["DataDictionary"]]

    ## currently we can only look at DataField  
    if(!all(names(dd) == "DataField")) warning("data specifications other than DataField are not yet implemented")
  
    ## check columns
    nc <- as.numeric(XML::xmlAttrs(dd)["numberOfFields"])
    if(!is.na(nc)) stopifnot(nc == length(dd))

    ## set up data frame (only numeric variables)
    mf <- as.data.frame(rep(list(1), nc))[0,]
    names(mf) <- XML::xmlSApply(dd, function(x) XML::xmlAttrs(x)["name"])

    ## modify class if necessary
    for(i in 1:nc) {
      optype <- XML::xmlAttrs(dd[[i]])["optype"]
      switch(optype,
        "categorical" = {
           mf[[i]] <- factor(integer(0),
	     levels = XML::xmlSApply(dd[[i]], function(x) gsub("&amp;", "&", XML::xmlAttrs(x)["value"], fixed = TRUE)))
        },
        "ordinal" = {
           mf[[i]] <- factor(integer(0), ordered = TRUE,
	     levels = XML::xmlSApply(dd[[i]], function(x) gsub("&amp;", "&", XML::xmlAttrs(x)["value"], fixed = TRUE)))
        },
        "continuous" = {
          dataType <- XML::xmlAttrs(dd[[i]])["dataType"]
          if(dataType == "integer") mf[[i]] <- integer(0)
        }
      )
    }
    
    return(mf)
  }
  mf <- extract_empty_model_frame(obj)
  mf_names <- names(mf)
  mf_levels <- lapply(mf, levels)

  ## parse MiningSchema
  extract_terms <- function(x) {

    ## extract MiningSchema
    ms <- x[["TreeModel"]]
    stopifnot("MiningSchema" %in% names(ms))
    ms <- ms[["MiningSchema"]]
    
    ## currently we can only look at MiningField  
    if(!all(names(ms) == "MiningField")) warning("MiningField not yet implemented")
    
    ## extract variable info
    vars <- t(XML::xmlSApply(ms, XML::xmlAttrs))
    if(sum(vars[,2] == "predicted") > 1) stop("multivariate responses not yet implemented")
    if(!all(vars[,2] %in% c("predicted", "active", "supplementary"))) warning("not yet implemented")
    
    ## set up formula
    ff <- as.formula(paste(vars[vars[,2] == "predicted",1], "~",
      paste(vars[vars[,2] != "predicted",1], collapse = " + ")))

    return(terms(ff))
  }
  trms <- extract_terms(obj)
  
  ## parse TreeModel
  tm <- obj[["TreeModel"]]
  tm_info <- c(XML::xmlAttrs(tm), h_info)

  ## check response
  stopifnot(tm_info["functionName"] %in% c("classification", "regression"))
  mf_response <- mf[[deparse(attr(trms, "variables")[[2L]])]]
  if(tm_info["functionName"] == "classification") stopifnot(inherits(mf_response, "factor"))
  if(tm_info["functionName"] == "regression") stopifnot(is.numeric(mf_response))
  
  ## convenience functions for parsing nodes
  is_terminal <- function(xnode) !("Node" %in% names(xnode))
  is_root <- function(xnode) "True" %in% names(xnode)
  n_kids <- function(xnode) sum("Node" == names(xnode))
  n_obs <- function(xnode) as.numeric(XML::xmlAttrs(xnode)["recordCount"])
  has_surrogates <- function(x) {
    ns <- sum(c("SimplePredicate", "SimpleSetPredicate", "CompoundPredicate") %in% names(x))
    if(ns != 1) stop("invalid PMML")
    if("CompoundPredicate" %in% names(x)) {
      if(identical(as.vector(XML::xmlAttrs(x[["CompoundPredicate"]])["booleanOperator"]), "surrogate")) return(TRUE)
        else return(FALSE)
    } else {
      return(FALSE)
    }
  }
  has_single_splits <- function(x) {
    wi <- which(names(x) %in% c("SimplePredicate", "SimpleSetPredicate", "CompoundPredicate"))
    sapply(wi, function(i) {
      if(names(x)[i] %in% c("SimplePredicate", "SimpleSetPredicate")) return(TRUE)
      if(identical(as.vector(XML::xmlAttrs(x[[i]])["booleanOperator"]), "or")) return(TRUE)
      stop("CompoundPredicate not yet implemented")
    })
  }
  n_splits <- function(xnode) {
    wi <- which("Node" == names(xnode))
    rval <- unique(sapply(wi, function(i) {
      xnodei <- if(has_surrogates(xnode[[i]])) xnode[[i]][["CompoundPredicate"]] else xnode[[i]]
      rval <- has_single_splits(xnodei)
      if(!all(rval)) stop("invalid PMML")
      sum(rval)
    }))
    if(length(rval) > 1) stop("invalid PMML")
    return(rval)
  }
  kid_ids <- function(xnode) {
    wi <- which("Node" == names(xnode))
    rval <- sapply(wi, function(j) {
      as.vector(XML::xmlAttrs(xnode[[j]])["id"])
    })  
  }
  get_pred <- function(xnode) {
    pred <- as.vector(XML::xmlAttrs(xnode)["score"])
    if(is.na(pred)) return(NULL)
    if(is.numeric(mf_response)) as.numeric(pred)
      else factor(pred, levels = levels(mf_response))
  }
  get_dist <- function(xnode) {
    wi <- which("ScoreDistribution" == names(xnode))
    if(length(wi) < 1) return(NULL)
    rval <- sapply(wi, function(i) as.numeric(XML::xmlAttrs(xnode[[i]])["recordCount"]))
    names(rval) <- sapply(wi, function(i) XML::xmlAttrs(xnode[[i]])["value"])
    if(inherits(mf_response, "factor")) rval <- rval[levels(mf_response)]
    return(rval)
  }
  get_error <- function(xnode) {
    if(tm_info["functionName"] != "classification") return(NULL)
    tab <- get_dist(xnode)
    if(is.null(tab)) return(NULL)
    c("%" = sum(100 * prop.table(tab)[names(tab) != get_pred(xnode)]))
  }
  get_extension <- function(xnode) {
    if(!("Extension" %in% names(xnode))) return(NULL)
    if(length(xnode[["Extension"]]) > 1) warning("currently only one Extension allowed")
    rval <- XML::xmlApply(xnode[["Extension"]][[1]], XML::xmlAttrs)
    names(rval) <- NULL
    rval <- unlist(rval)
    to_numeric <- function(x) {
      y <- suppressWarnings(as.numeric(x))
      if(!is.null(y) && !is.na(y)) y else x
    }
    sapply(rval, to_numeric)
  }
  node_info <- function(xnode) list(prediction = get_pred(xnode), n = n_obs(xnode),
    error = get_error(xnode), distribution = get_dist(xnode), extension = get_extension(xnode))
  get_split_prob <- function(xnode) {
    rval <- rep(0, n_kids(xnode))
    wi <- XML::xmlAttrs(xnode)["defaultChild"]
    if(is.na(wi)) rval <- NULL
      else rval[which(kid_ids(xnode) == wi)] <- 1
    return(rval)
  }
  get_split <- function(xnode, i, surrogates) {
    wi <- which("Node" == names(xnode))
    rval <- sapply(wi, function(j) {
      nj <- if(surrogates) xnode[[j]][["CompoundPredicate"]] else xnode[[j]]
      if(any(c("SimplePredicate", "SimpleSetPredicate") %in% names(nj))) {
        wii <- which(names(nj) %in% c("SimplePredicate", "SimpleSetPredicate"))[i]
        c("predicateType" = as.vector(names(nj)[wii]), XML::xmlAttrs(nj[[wii]]))
      } else {
        wii <- which(names(nj) == "CompoundPredicate")[i]
	nj <- nj[[wii]]
	if(!identical(as.vector(XML::xmlAttrs(nj)["booleanOperator"]), "or")) stop("not yet implemented")
	if(any(names(nj) %in% c("SimpleSetPredicate", "CompoundPredicate"))) stop("not yet implemented")
	rvali <- sapply(which(names(nj) == "SimplePredicate"), function(j)
	  c("predicateType" = as.vector(names(nj)[j]), XML::xmlAttrs(nj[[j]])))
	if(is.null(dim(rvali))) rvali <- matrix(rvali, ncol = 1)
	stopifnot(length(unique(rvali["predicateType",])) == 1)
	stopifnot(length(unique(rvali["field",])) == 1)
	stopifnot(all(rvali["operator",] == "equal"))
	c("predicateType" = "simpleSetPredicate",
	  "field" = rvali["field", 1],
	  "booleanOperator" = "isIn")
      }
    })
    stopifnot(length(unique(rval["predicateType",])) == 1)
    stopifnot(length(unique(rval["field",])) == 1)    
    if(rval["predicateType", 1] == "SimplePredicate") {
      stopifnot(length(unique(rval["value",])) == 1)
      if(ncol(rval) != 2) stop("not yet implemented")
      if(!(identical(as.vector(sort(rval["operator",])), c("greaterThan", "lessOrEqual")) |
           identical(as.vector(sort(rval["operator",])), c("greaterOrEqual", "lessThan")))
      ) stop("not yet implemented")
      partysplit(
        varid = which(rval["field", 1] == mf_names),
	breaks = as.numeric(rval["value", 1]),
	index = if(substr(rval["operator", 1], 1, 1) != "l") 2:1 else NULL,
	right = "lessOrEqual" %in% rval["operator",],
	prob = if(i == 1) get_split_prob(xnode) else NULL
      )      
    } else {
      varid <- which(rval["field", 1] == mf_names)
      lev <- mf_levels[[varid]]
      stopifnot(length(lev) > 1)
      idx <- rep(NA, length(lev))
      lab <- lapply(wi, function(j) {
        nj <- if(surrogates) xnode[[j]][["CompoundPredicate"]] else xnode[[j]]
        if(any(names(nj) %in% c("SimplePredicate", "SimpleSetPredicate"))) {
          wii <- which(names(nj) %in% c("SimplePredicate", "SimpleSetPredicate"))[i]      
          ar <- nj[[wii]][["Array"]]
	  stopifnot(XML::xmlAttrs(ar)["type"] == "string")
	  rv <- XML::xmlValue(ar)
	  rv <- gsub("&quot;", "\"", rv, fixed = TRUE)
	  rv <- if(substr(rv, 1, 1) == "\"" & substr(rv, nchar(rv), nchar(rv)) == "\"") {
	    strsplit(substr(rv, 2, nchar(rv) - 1), "\" \"")[[1]]
	  } else {
	    strsplit(rv, " ")[[1]]
	  }
	  stopifnot(length(rv) == as.numeric(XML::xmlAttrs(ar)["n"]))
	  return(rv)
	} else {
          wii <- which(names(nj) == "CompoundPredicate")[i]	
	  as.vector(XML::xmlSApply(nj[[wii]], function(z) XML::xmlAttrs(z)["value"]))
	}
      })
      for(j in 1:ncol(rval)) {
        if(rval["booleanOperator",j] == "isIn") idx[which(lev %in% lab[[j]])] <- j
	  else idx[which(!(lev %in% lab[[j]]))] <- j
      }
      stopifnot(all(na.omit(idx) > 0))
      if(min(idx, na.rm = TRUE) != 1) stop(sprintf("variable levels (%s) and split labels (%s)",
        paste(lev, collapse = ", "), paste(sapply(lab, paste, collapse = ", "), collapse = " | ")))
      
      partysplit(
        varid = varid,
	breaks = NULL,
	index = as.integer(idx),
	prob = if(i == 1) get_split_prob(xnode) else NULL
      )
    }
  }
  
  ## function for setting up nodes
  ## (using global index ii)
  pmml_node <- function(xnode) {
    ii <<- ii + 1
    if(is_terminal(xnode)) return(partynode(as.integer(ii),
      info = node_info(xnode)
    ))
    wi <- which("Node" == names(xnode))
    ns <- n_splits(xnode)    
    nd <- partynode(as.integer(ii),
      split = get_split(xnode, 1, has_surrogates(xnode[[wi[1]]])),
      kids = lapply(wi, function(j) pmml_node(xnode[[j]])),
      surrogates = if(ns < 2) NULL else lapply(2:ns, function(j) get_split(xnode, j, TRUE)),
      info = node_info(xnode)
    )
    nd
  }
  
  ## set up node
  ii <- 0
  if(is_root(tm[["Node"]])) nd <- pmml_node(tm[["Node"]]) else stop("root node not declared, invalid PMML?")

  ## set up party
  ## FIXME: extend info slot?
  pt <- party(node = nd, data = mf, fitted = NULL, terms = trms, names = NULL, info = tm_info)
  class(pt) <- c("simpleparty", class(pt))

  return(pt)
}
