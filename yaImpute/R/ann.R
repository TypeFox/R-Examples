"ann" <- function(ref, target, k=1, eps=0.0, tree.type="kd", search.type="standard", bucket.size=1, split.rule="sl_midpt", shrink.rule="simple", verbose=TRUE, ...){

  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }
  
  if(missing(ref)){stop("error: ref must be specified")}
  if(missing(target)){stop("error: target must be specified")}

  if(! is.matrix(ref)){stop("error: ref must be matrix")}
  if(! is.matrix(target)){stop("error: target must be matrix")}
  if(k <= 0){stop("error: k must be int > 0")}
  if(eps < 0.0){stop("error: eps must be int >= 0.0")}
  if(!tree.type %in% c("kd", "bd", "brute")){stop("error: tree.type must be 'kd', 'bd', or 'brute'")}
  if(!search.type %in% c("standard", "priority")){stop("error: tree.type must be 'standard' or 'priority'")}

  pri <- FALSE
  if(search.type == "priority") pri <- TRUE

  if(nrow(ref) == 0 || nrow(target) == 0){stop("error: nrow(ref) and nrow(target) must be > 0")}
  if(ncol(ref) != ncol(target)){stop("error: ncol(ref) must equal ncol(target)")}
  if(k > nrow(ref)){stop("error: k must be <= nrow(ref)")}
  if(bucket.size <=0){stop("error: bucket.size must be > 0")}

  split <- list("standard"=0, "midpt"=1, "fair"=2, "sl_midpt"=3, "sl_fair"=4, "suggest"=5)
  if(!split.rule %in% names(split)){stop("error: ",split.rule ," is not a valid split rule, choices are 'standard', 'midpt', 'fair', 'sl_midptm', 'sl_fair', and 'suggest'")}
  split.rule <- split[[split.rule]]
  
  shrink <- list("none"=0, "simple"=1, "centroid"=2, "suggest"=3)
  if(!shrink.rule %in% names(shrink)){stop("error: ",shrink.rule ," is not a valid shrink rule for the bd-tree, choices are 'none', 'simple', 'centroid', and 'suggest'")}
  shrink.rule <- shrink[[shrink.rule]]
  
  storage.mode(ref) <- storage.mode(target) <- "double"

  if(tree.type == "bd" && any(duplicated(ref))){stop("error: duplicate pattern found in the 'ref' matrix; therefore, bd-tree cannot be used")}
  
  args <- list("ref"=ref, "target"=target, "k"=as.integer(k), "eps"=as.double(eps), "tree.type"=tree.type,
               "priority"=as.integer(pri), "bucket.size"=as.integer(bucket.size),
               "split.rule"=as.integer(split.rule), "shrink.rule"=as.integer(shrink.rule), "verbose"=as.integer(verbose))
  
  out <- .Call("ann", args)
  out$k <- k
  out$tree.type <- tree.type

  if(tree.type %in% c("kd", "bd")){
    out$eps <- eps
    out$search.type <- search.type
    out$bucket.size <- bucket.size
    out$split.rule <- split.rule
    if(tree.type == "bd")
      out$shrink.rule <- shrink.rule
  }
  
  class(out) <- "ann"
  out
}
