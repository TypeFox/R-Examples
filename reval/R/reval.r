#' This package provides the function \code{evalmany} to simplify scenario 
#' testing and sensitivity analysis with R.
#'
#' @name reval-package
#' @docType package
#' @import foreach parallel doParallel
NULL


eval_set = function(fun, pargs, default.args, pkg){
  recycle.grid = function(...){
    dotList = list(...)
    max.length = max(sapply(dotList, length))
    lapply(dotList, rep, length = max.length)
  }
  # setup inputs
  setpargs = do.call(recycle.grid, pargs) # recycle values for sets  
  for(n in names(setpargs))
    if(is.null(names(setpargs[[n]])))
      names(setpargs[[n]]) = paste(setpargs[[n]])
  numsets = length(setpargs[[1]])
  numpargs = length(setpargs)
  # setup outputs
  n = vector("list", length = numsets)
  for(i in seq(numsets))
    n[[i]] = vector("list", length = numpargs)
  setnames = rep(NA, numsets)
  for(j in seq(numsets)){
    for(i in seq(numpargs))
      n[[j]][[i]] = paste(names(setpargs)[i], names(setpargs[[i]])[j], 
        sep = " = ")
    setnames[j] = paste(n[[j]], collapse = " ; ")
    }
  # prepare argument sets
  argsets = vector("list", length = numsets)
  names(argsets) = setnames  
  for(i in seq(numsets)){
   args = default.args
     for(p in names(setpargs))
       args[[p]] = setpargs[[p]][[i]]
   argsets[[setnames[i]]] = args
  }  
  # evaluate
  an = NULL # stupidity for CRAN
  as = NULL
  ret = foreach(as = argsets, an = setnames, .errorhandling = "pass", 
    .inorder = FALSE, .combine = 'c', .packages = pkg) %dopar% {
    tryCatch({
      structure(list(do.call(fun, as)), names = an)
    }, error = function(cond){
      structure(list(cond), names = an)
    }, warning = function(cond){
      structure(list(cond), names = an)
    })
  }
  # error handling
  for(n in names(ret))
    if(any(c("error", "warning") %in% class(ret[[n]]))){
      message(ret[[n]], " when evaluating ", n)
      if("error" %in% class(ret[[n]]))
        ret[[n]] = NULL
    }
  return(ret)
}

eval_permute = function(fun, pargs, default.args, size, pkg){
  arg.grid = expand.grid(pargs)
  if(size < 1)
    size = nrow(arg.grid)  
  permargs = as.list(arg.grid[sample(seq(nrow(arg.grid)), size = size),])
  attributes(permargs) = NULL 
  names(permargs) = names(arg.grid)
  eval_set(fun, permargs, default.args, pkg)
}

eval_ofat = function(fun, pargs, default.args, pkg){
  base.args = lapply(pargs, function(x) x[1])
  ret = list()
  for(p in names(pargs)){
    thisargs = base.args
    thisargs[[p]] = pargs[[p]]
    ret = append(ret, eval_set(fun, thisargs, default.args, pkg))
  }
  return(ret)
}

collate_set = function(l, collate.fun, collate.id, prepend){
  id.name = paste0(prepend, "id")
  res = NULL
  for(n in names(l))
    if(!is.null(l[[n]])){
      nres = as.data.frame(collate.fun(l[[n]]))
      nres[id.name] = n
      res = rbind(res, nres)
    }
  if(collate.id == "single")
    res
  else
    parse_set(res, id.name, prepend)
}

parse_set = function(df, n, prepend){
  x = df[[n]]
  l = NULL
  variables = strsplit(x, " ; ") 
  for(i in seq(length(variables))){
    value = strsplit(variables[[i]], " = ")
    v = sapply(value, function(x) x[2])
    names(v) = paste0(prepend, sapply(value, function(x) x[1]))
    l = rbind(l, v)
  }
  rownames(l) = NULL
  out = cbind(df, l)
  out[[n]] = NULL
  out
}

#' Repeated evaluations
#'
#' Evaluate a function repeatedly across argument sets or permutations.
#'
#' @param fun The function to be evaluated.
#' @param ... Arguments to be varied when evaluating \code{fun}, where each 
#'   argument in '\code{...}' is a (named) vector or list of values. Lists of
#'   multi-value objects (e.g. data.frames) should be named explicitly and may 
#'   otherwise produce unexpected or incorrect names.
#' @param method The sensitivity analysis method to be used. Can be either 
#'   one-factor-at-a-time ("ofat") evaluation, evaluation of parameter sets
#'   ("set"), or (sampled) permutations of parameter sets ("permute"). When
#'   \code{method = "ofat"}, the first element of each argument in '\code{...}'
#'   is assumed to be the "default" value of that argument.
#' @param sample If \code{method = "permute"}, the number of parameter 
#'   permutations to evaluate (sampling without replacement). 
#'   If \code{sample < 1} (the default) then all possible permutations are 
#'   evaluated.
#' @param default.args The default values of arguments passed to \code{fun}.
#' @param collate Whether to collate the results or not. If TRUE, output 
#'   elements will be coerced into data.frames using \code{as.data.frame}.
#'   Otherwise, the raw outputs will be returned as a named list.
#' @param collate.id If \code{collate = TRUE}, the method used to 
#'   store the evaluation identifiers. If \code{collate.id = "single"}, a 
#'   single column named 'id' is used. If \code{collate.id = "multi"}, 
#'   one column is created for each argument in '\code{...}', e.g. 
#'   'arg1', 'arg2', etc. 
#' @param collate.prepend A character string prepended to the identifier 
#'   column. If \code{collate.id = "single"}, the identifier column will be
#'   named \code{<collate.prepend>id}. If \code{collate.id = "multi"}, 
#'   identifier columns will be named as \code{<collate.prepend><arg>} where
#'   \code{arg} is an element of \code{...}.
#' @param collate.fun If \code{collate = TRUE}, an optional function 
#'   for reshaping the output of each evaluation prior to coercing and 
#'   collating the outputs. 
#' @param clusters Number of clusters to use for parallel processing. Default
#'   is 1 (serial computation).
#' @param packages For parallel processing. Character vector of packages that 
#'   'fun' depends on.
#'
#' @return If \code{collate = TRUE}, a data.frame. Otherwise, a named list.
#'
#' @examples
#' myfun = function(n, mean=0, sd = 1){ 
#'   x = rnorm(n, mean, sd) 
#'   data.frame(sample.mean = mean(x), sample.sd = sd(x))
#' }
#' evalmany(myfun, mean = c(5, 9), sd = c(2, 3), default.args = list(n = 1e6))
#' evalmany(myfun, mean = seq(20), sd = seq(1, 4, by = 0.1), 
#'   default.args = list(n = 1e6), method = "permute", sample = 10,
#'   collate.id = "multi", collate.prepend = "arg.")
#'
#' # vector recycling
#' evalmany(myfun, mean = c(0, 3, 5), sd = c(1, 10), 
#'   default.args = list(n = 1e6), method = "set", collate.id = "multi")
#' 
#' # Parallel processing
#' evalmany(myfun, mean = seq(0, 50, by = 10), sd = seq(1, 10, by = 1.5), 
#'   default.args = list(n = 1e5), method = "permute", collate.id = "multi",  
#'   clusters = 2)
#'
#' \dontrun{
#' # Complex objects 
#' formulas = list(y ~ 1, y ~ x, y ~ x + z)   
#' datasets = list(
#'   A = data.frame(x = seq(0, 99), y = seq(0, 99) + rnorm(100)),
#'   B = data.frame(x = seq(0, 99), y = seq(0, 99) + rnorm(100, mean = 5)),
#'   C = data.frame(x = seq(0, 99), y = seq(0, 99) + rlnorm(100, meanlog = 1),
#'     z = seq(0, 99) - rlnorm(100, meanlog = -1))
#' )
#' # raw output
#' evalmany(lm, formula = formulas, data = datasets, method = "set", 
#'   collate = FALSE)
#' # data extraction and error handling
#' evalmany(lm, formula = formulas, data = datasets, method = "permute",
#'   collate.id = "multi", collate.fun = function(x) 
#'     data.frame(param = names(x$coefficients), value = x$coefficients, 
#'     row.names=NULL))
#' }
#'
#' @export
evalmany = function(fun, ..., method = c("ofat", "permute", "set"), 
 sample = 0L, default.args = list(), collate = TRUE, 
 collate.id = c("single", "multi"), collate.prepend = "", 
 collate.fun = identity, clusters = 1L, packages = NULL){
  # argument checks
  if(any(names(list(...)) == ""))
    stop("All arguments to 'fun' must be named")
  method = match.arg(method, c("ofat", "permute", "set"))
  collate.id = match.arg(collate.id, c("single", "multi"))
  pargs = list(...)
  # name sanitation
  if(any(sapply(pargs, function(x) any(grepl(";", names(x))))))
    stop('argument names must not contain ";".')
  # prepare for parallel processing
  if(clusters > 1){
    cl <- makeCluster(clusters)
    on.exit(stopCluster(cl))
    registerDoParallel(cl)
  }
  else
    registerDoSEQ()
  # process
  if(method == "ofat")
    res = eval_ofat(fun, pargs, default.args, packages)
  else if(method == "set")
    res = eval_set(fun, pargs, default.args, packages)
  else
    res = eval_permute(fun, pargs, default.args, size = as.integer(sample), 
      packages)
  # format outputs
  if(collate == FALSE)
    res
  else
    collate_set(res, collate.fun, collate.id, collate.prepend)
}
