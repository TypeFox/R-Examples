### Copied from R-3.0.1, to support getKnownS3generics.
getKnownS3generics <- function(){
  c(names(.knownS3Generics), get_internal_S3_generics())
}

### Copied from R-3.0.1, to support getKnownS3generics.
get_internal_S3_generics <- function(primitive = TRUE){
  out <- c("[", "[[", "$", "[<-", "[[<-", "$<-", "as.vector", 
           "unlist", get_S3_primitive_generics())
  if (!primitive) 
    out <- out[!vapply(out, is_primitive_in_base, NA)]
  out
}

### Copied from R-3.0.1, to support getKnownS3generics.
is_primitive_in_base <- function(fname){
  is.primitive(get(fname, envir = baseenv(), inherits = FALSE))
}

### Copied from R-3.0.1, to support getKnownS3generics.
get_S3_primitive_generics <- function(include_group_generics = TRUE){
  if (include_group_generics) 
    c(base::.S3PrimitiveGenerics, "abs", "sign", "sqrt", 
      "floor", "ceiling", "trunc", "round", "signif", "exp", 
      "log", "expm1", "log1p", "cos", "sin", "tan", "acos", 
      "asin", "atan", "cosh", "sinh", "tanh", "acosh", 
      "asinh", "atanh", "lgamma", "gamma", "digamma", "trigamma", 
      "cumsum", "cumprod", "cummax", "cummin", "+", "-", 
      "*", "/", "^", "%%", "%/%", "&", "|", "!", "==", 
      "!=", "<", "<=", ">=", ">", "all", "any", "sum", 
      "prod", "max", "min", "range", "Arg", "Conj", "Im", 
      "Mod", "Re")
  else base::.S3PrimitiveGenerics
}

### Copied from R-3.0.1, to support findGeneric.
findGeneric <- function(fname, envir){
  if (!exists(fname, mode = "function", envir = envir)) 
    return("")
  f <- get(fname, mode = "function", envir = envir)
  if (.isMethodsDispatchOn() && methods::is(f, "genericFunction")) {
    fMethsEnv <- methods::getMethodsForDispatch(f)
    r <- lapply(grep("^ANY\\b", ls(envir = fMethsEnv), value = TRUE), 
                get, envir = fMethsEnv)
    if (any(ddm <- unlist(lapply(r, class)) == "derivedDefaultMethod")) 
      f <- r[ddm][[1]]@.Data
    else warning(gettextf("'%s' is a formal generic function; S3 methods will not likely be found", 
                          fname), domain = NA)
  }
  isUMEbrace <- function(e) {
    for (ee in as.list(e[-1L])) if (nzchar(res <- isUME(ee))) 
      return(res)
    ""
  }
  isUMEif <- function(e) {
    if (length(e) == 3L) 
      isUME(e[[3L]])
    else {
      if (nzchar(res <- isUME(e[[3L]]))) 
        res
      else if (nzchar(res <- isUME(e[[4L]]))) 
        res
      else ""
    }
  }
  isUME <- function(e) {
    if (is.call(e) && (is.name(e[[1L]]) || is.character(e[[1L]]))) {
      switch(as.character(e[[1L]]), UseMethod = as.character(e[[2L]]), 
             `{` = isUMEbrace(e), `if` = isUMEif(e), "")
    }
    else ""
  }
  isUME(body(f))
}

### Copied from R-3.0.1, to support fixPackageFileNames.
fixPackageFileNames <- function(list){
  list <- as.character(list)
  if (length(list) == 0L) 
    return(list)
  list0 <- gsub("[[:cntrl:]\"*/:<>?\\|]", "_", list)
  wrong <- grep("^(con|prn|aux|clock\\$|nul|lpt[1-3]|com[1-4])(\\..*|)$", 
                list0)
  if (length(wrong)) 
    list0[wrong] <- paste0("zz", list0[wrong])
  ok <- grepl("^[[:alnum:]]", list0)
  if (any(!ok)) 
    list0[!ok] <- paste0("z", list0[!ok])
  list1 <- tolower(list0)
  list2 <- make.unique(list1, sep = "_")
  changed <- (list2 != list1)
  list0[changed] <- list2[changed]
  list0
}

