# internals.r: Various utility functions which are (usually) not called from
# outside the package


### Functions to create SQL statements and to retreive record pairs from
### the database

# construct part of where-clause which represents blocking restrictions
# make something like 'list(1:2, 3)' to something like
# '(t1.field1=t2.field1 and t1.field2=t2.field2) or (t1.field3=t2.field3)'
blockfldfun <- function(blockfld, phoneticFld, phoneticFun, coln)
{
  blockElemFun <- function(fldIndex)
  {
    if (fldIndex %in% phoneticFld)
      return(sprintf("%1$s(t1.'%2$s')=%1$s(t2.'%2$s')", phoneticFun, coln[fldIndex]))
    else
      return(sprintf("t1.'%1$s'=t2.'%1$s'", coln[fldIndex]))
  }

 paste("(", paste(sapply(blockfld, function(blockvec)
                  paste(sapply(blockvec, blockElemFun),
                        collapse=" and ")),
                  collapse=") or ("), ")", sep="")
}

#' Create SQL statement
#'
#' Creates SQL statememt to retreive comparison patterns, respecting
#' parameters such as blocking definition and exclusion of fields.
#'
#' @value A list with components "select_list", "from_clause", "where_clause"
#' representing the corresponding parts of the query without the keywords
#' 'SELECT', 'FROM' and 'WHERE'.
getSQLStatement <- function(data1, data2 = data1, con, type, blockFld,
  excludeFld, strcmpFld, strcmpFun, phoneticFld, phoneticFun)
{
  # constructs select for a single column, to be used by lapply 
  # (see below)
  selectListElem <- function(fldIndex, coln, excludeFld, strcmpFld, strcmpFun,
                            phoneticFld, phoneticFun)
  {
    # nothing if field is excluded
    if (fldIndex %in% excludeFld)
      return(character(0))
      
    # enclose fields in phonetic function if desired
    if (fldIndex %in% phoneticFld)
    {
      fld1 <- sprintf("%s(t1.'%s')", phoneticFun, coln[fldIndex])
      fld2 <- sprintf("%s(t2.'%s')", phoneticFun, coln[fldIndex])
    } else
    {
      fld1 <- sprintf("t1.'%s'", coln[fldIndex])
      fld2 <- sprintf("t2.'%s'", coln[fldIndex])
    }


    # something like 'jarowinkler(t1.fname, t2.fname) as fname'
    if (fldIndex %in% strcmpFld)
      return(sprintf("%s(%s, %s) as '%s'", strcmpFun, fld1, fld2, coln[fldIndex]))


    # direct comparison: something like 't1.fname=t2.fname as fname'      
    return(sprintf("%s=%s as '%s'", fld1, fld2, coln[fldIndex]))
  }
  coln <- make.db.names(con, colnames(data1), keywords = SQLKeywords(con))

  selectlist_id <- "t1.row_names as id1, t2.row_names as id2"
  # use unlist to delete NULLs from list
  selectlist <- paste(unlist(lapply(1:length(coln), selectListElem,
    coln, excludeFld, strcmpFld, strcmpFun,
    phoneticFld, phoneticFun)), collapse = ", ")
  selectlist <- paste(selectlist, "t1.identity=t2.identity as is_match", sep=",")
  fromclause <- switch(type, deduplication = "data t1, data t2",
                                      linkage = "data1 t1, data2 t2")
  whereclause <- switch(type, deduplication = "t1.row_names < t2.row_names",
                                      linkage = "1")
  if (length(blockFld)>0)
  {
   whereclause <- sprintf("%s and (%s)", whereclause, blockfldfun(blockFld,
    phoneticFld, phoneticFun, coln))
  }
  return(list(select_list = paste(selectlist_id, selectlist, sep=", "),
                from_clause = fromclause, where_clause = whereclause))
}


#' Begin generation of data pairs
#'
#' An SQL statement representing the generation of data pairs, including
#' the configuration of blocking fields, phonetics etc. is constructed and
#' send to SQLite.
setGeneric(
  name = "begin",
  def = function(x, ...) standardGeneric("begin")
)
      
setMethod(
  f = "begin",
  signature = "RLBigData",
  definition = function(x, ...)
  {
    sql <- getSQLStatement(x)  
    query <- sprintf("select %s from %s where %s", sql$select_list, 
      sql$from_clause, sql$where_clause)
    dbSendQuery(x@con, query) # can be retreived via dbListResults(x@con)[[1]]
    invisible(x)
  }
)

# retreive next n pairs
setGeneric(
  name = "nextPairs",
  def = function(x, n=10000, ...) standardGeneric("nextPairs")
)

setMethod(
  f = "nextPairs",
  signature = "RLBigData",
  definition = function(x, n=10000, ...)
  {
    res <- dbListResults(x@con)[[1]]
    result <- fetch(res, n)
    # Spalten, die nur NA enthalten, werden als character ausgegeben, deshalb
    # Umwandlung nicht-numerischer Spalten in numeric
    if (nrow(result) > 0) # wichtig, weil sonst Fehler bei Zugriff auf Spalten auftritt
    {
      for (i in 1:ncol(result))
      {
        if (!is.numeric(result[,i]))
          result[,i] <- as.numeric(result[,i])
      }
    }
    # column names may have changed due to SQL conform conversion, reset them
    colnames(result)[-c(1,2,ncol(result))] <- getColumnNames(x)
    result
  }
)

# Clean up after retreiving data pairs (close result set in database)
setGeneric(
  name = "clear",
  def = function(x, ...) standardGeneric("clear")
)

setMethod(
  f = "clear",
  signature = "RLBigData",
  definition = function(x, ...) dbClearResult(dbListResults(x@con)[[1]])
)


### Functions neccessary to load extensions into SQLite database

# Function body taken from init_extension, package RSQLite.extfuns
init_sqlite_extensions <- function(db)
{
    ans <- FALSE
    if (.allows_extensions(db)) {
        res <- dbGetQuery(db, sprintf("SELECT load_extension('%s')",
                                      .lib_path()))
        ans <- all(dim(res) == c(1, 1))
    } else {
        stop("loadable extensions are not enabled for this db connection")
    }
    ans
}

# taken from RSQLite.extfuns
.allows_extensions <- function(db)
{
    v <- db@loadable.extensions
    isTRUE(v)
}


# taken from RSQLite.extfuns
.lib_path <- function()
{
    ## this is a bit of a trick, but the NAMESPACE code
    ## puts .packageName in the package environment and this
    ## seems slightly better than hard-coding.
    ##
    ## This also relies on the DLL being loaded even though only SQLite
    ## actually needs to load the library.  It does not appear that
    ## loading it causes any harm and it makes finding the path easy
    ## (don't have to worry about arch issues).
    getLoadedDLLs()[[.packageName]][["path"]]
}



# get count of each distinct comparison pattern
# Fuzzy values above cutoff are converted to 1, below cutoff to 0
# NAs are converted to 0
setGeneric(
  name = "getPatternCounts",
  def = function(x, n=10000, cutoff=1, withProgressBar = (sink.number()==0))
    standardGeneric("getPatternCounts")
)

setMethod(
  f = "getPatternCounts",
  signature = "RLBigData",
  definition = function(x, n=10000, cutoff=1, withProgressBar = (sink.number()==0))
  {
    if (withProgressBar)
    {
      pgb <- txtProgressBar(max=nrow(x@pairs))
    }

    patternCounts <- 0L
    ffrowapply(
      {
        slice <- as.ram(x@pairs[i1:i2, 3:(ncol(x@pairs)-1), drop=FALSE])
        slice[is.na(slice)] <- 0
        slice[slice < cutoff] <- 0
        slice[slice >= cutoff] <- 1
        patternCounts <- patternCounts + countpattern(slice)
        if (withProgressBar) setTxtProgressBar(pgb, i2)
      }, X = x@pairs)
      patternCounts
  }
)

# get number of matches
setGeneric(
  name = "getMatchCount",
  def = function(object) standardGeneric("getMatchCount")
)

setMethod(
  f = "getMatchCount",
  signature = "RLBigData",
  definition = function(object)
  {
    sum(object@pairs$is_match, na.rm=TRUE)
  }
)

# get number of non-matches
setGeneric(
  name = "getNonMatchCount",
  def = function(object) standardGeneric("getNonMatchCount")
)

setMethod(
  f = "getNonMatchCount",
  signature = "RLBigData",
  definition = function(object)
  {
    sum(chunkify(function(x) x==0)(object@pairs$is_match), na.rm=TRUE)
  }
)

# Get the number of pairs with unknown matching status
setGeneric(
  name = "getNACount",
  def = function(object) standardGeneric("getNACount")
)

setMethod(
  f = "getNACount",
  signature = "RLBigData",
  definition = function(object)
  {
    sum(chunkify(is.na)(object@pairs$is_match))
  }
)


setGeneric(
  name = "getColumnNames",
  def = function(object, withExcluded = FALSE) standardGeneric("getColumnNames")
)

setMethod(
  f = "getColumnNames",
  signature = "RLBigDataDedup",
  definition = function(object, withExcluded = FALSE)
  {
    if (withExcluded || length(object@excludeFld)==0) colnames(object@data)
    else colnames(object@data)[-object@excludeFld]
  }
)

setMethod(
  f = "getColumnNames",
  signature = "RLBigDataLinkage",
  definition = function(object, withExcluded = FALSE)
  {
    if (withExcluded || length(object@excludeFld)==0) colnames(object@data1)
    else colnames(object@data1)[-object@excludeFld]
  }
)

### Various other utility functions

# utility function to generate all unordered pairs of x[1]..x[n]
# if x is a vector, or 1..x
unorderedPairs <- function (x)
{
    if (length(x)==1)
    {
      if (!is.numeric(x) || x < 2)
        stop("x must be a vector or a number >= 2")
        return (array(unlist(lapply(1:(x-1),
          function (k) rbind(k,(k+1):x))),dim=c(2,x*(x-1)/2)))
    }
    if (!is.vector(x))
      stop ("x must be a vector or a number >= 2")
    n=length(x)
    return (array(unlist(lapply(1:(n-1),
    function (k) rbind(x[k],x[(k+1):n]))),dim=c(2,n*(n-1)/2)))
}

isFALSE <- function(x) identical(x,FALSE)

deleteNULLs  <-  function(x)
    x[unlist(lapply(x, length) != 0)]

# interprets a scalar x as a set with 1 element (see also man page for sample)
resample <- function(x, size, ...)
     if(length(x) <= 1) { if(!missing(size) && size == 0) x[FALSE] else x
     } else sample(x, size, ...)


# modified function from package e1071, works also for data with only one column
countpattern <- function (x, matching = FALSE)
{
    nvar <- dim(x)[2]
    n <- dim(x)[1]
    b <- matrix(0, 2^nvar, nvar)
    for (i in 1:nvar) b[, nvar + 1 - i] <- rep(rep(c(0, 1), c(2^(i -
        1), 2^(i - 1))), 2^(nvar - i))
    namespat <- character(nrow(b))
    for (i in 1:nvar) namespat <- paste(namespat, b[, i], sep = "")
    xpat <- numeric(nrow(x))
    for (i in 1:nvar) xpat <- 2 * xpat + x[, i]
    xpat <- xpat + 1
    pat <- tabulate(xpat, nbins = 2^nvar)
    names(pat) <- namespat
    if (matching)
        return(list(pat = pat, matching = xpat))
    else return(pat)
}
