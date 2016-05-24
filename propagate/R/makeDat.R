makeDat <- function(expr) { 
  
  if (!is.expression(expr)) stop("input must be an expression!")
  
  ## get all variables names from workspace
  LS <- ls(envir = .GlobalEnv)
  
  ## get all variable names from expression
  exprNAMES <- all.vars(expr)  
  
  ## match both
  MATCH <- match(all.names(expr), LS)
  SEL <- unique(MATCH[!is.na(MATCH)])
  lsNAMES <- LS[SEL]
  
  ## check for presence of identical variable names
  ## in expression and workspace
  if (length(exprNAMES) != length(lsNAMES)) {
    DIFF <- setdiff(exprNAMES, lsNAMES)
    stop("object '", DIFF, "' in expresssion but not in workspace!")    
  }
  
  ## create data list and combine to dataframe
  LIST <- lapply(exprNAMES, function(x) get(x, envir = .GlobalEnv))
  DATA <- do.call("cbind", LIST)
  colnames(DATA) <- exprNAMES
  
  return(DATA)  
}
