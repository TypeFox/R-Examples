
# Check an input table
checkTbl= function(tbl, tblName, colNames, nameCol, emptyOK) {
  # Check type
  if (!is.data.frame(tbl))
    stop(paste0("'",tblName,"' must be a data frame"))
  # Handle empty tables
  if ((nrow(tbl) == 0) && (!emptyOK)) {
    stop(paste0("number of records in '",tblName,"' must not be zero"))
  } else if ((nrow(tbl) == 0) && emptyOK) {
    return(NULL)
  # Handle tables with some contents
  } else {
    for (n in colNames)
      tbl[,n]= as.character(tbl[,n])
    # Check column names
    if (!all(colNames %in% names(tbl)))
      stop(paste0("'",tblName,"' must have columns '",
        paste(colNames,collapse="', '"),"'"))
    # Check entries in name column
    if (!is.null(nameCol)) {
      if (!(nameCol %in% names(tbl)))
        stop(paste0("column '",nameCol," not present in '",tblName,"'"))
      if (any(is.na(tbl[,nameCol])))
        stop(paste0("NA values not allowed in column '",nameCol," of '",
          tblName,"'"))
      if (any(duplicated(tbl[,nameCol])))
        stop(paste0("duplicate names detected in column '",nameCol," of '",
          tblName,"'"))
      if (any(tbl[,nameCol] %in% rodeoConst$reservedNames))
        stop(paste0("column '",nameCol," of '",tblName,"' must not",
          " contain any of the reserved words '",
          paste(rodeoConst$reservedNames, collapse="', '"),"'"))
    }
    # Check whether names are legal identifiers
    bad= tbl[,nameCol][!grepl(pattern="^[a-zA-Z]+[a-zA-Z0-9_]*$", x=tbl[,nameCol])]
    if (length(bad) > 0)
      stop(paste0("the following name(s) from column '",nameCol,"' of '",tblName,
        "' are not legal identifiers: '",paste(bad,collapse="', '"),"'"))
  }
  return(NULL)
}

# Extract identifiers from a mathematical expression (given as a string)
extractIdentifiers= function(expr, dropDuplicates=TRUE) {
  tmp= gregexpr(pattern=rodeoConst$identifierPatterns$core,text=expr)[[1]]
  if (tmp[1] == -1) {
    return(character(0))
  } else {
    first= tmp
    final= tmp-1+attr(tmp,which="match.length",exact=TRUE)
    res= substring(expr, first, final)
    return(ifelse(dropDuplicates,unique(res),res))
  }
}

# Find undeclared identifiers in a mathematical expression (given as a string)
undeclared= function(expr, knownNames) {
  ident= extractIdentifiers(expr)
  ident[!(ident %in% knownNames)]
}

# Substitute identifiers in a mathematical expression (given as a string)
substituteIdentifiers= function(expr, sub, all=TRUE) {
  # Check vector of substitutes
  if (is.null(names(sub)) || any(names(sub) == ""))
    stop("missing element name(s) in vector of substitutes")
  if (any(duplicated(names(sub))))
    stop("duplicated element name(s) in vector of substitutes")
  if (any((names(sub) %in% sub) & (names(sub) != sub)))
    stop("bad vector of substitutes (the VALUE of an element must not be",
      "identical to the NAME of another element)")
  specialChar="\a"
  if (grepl(pattern=specialChar, x=expr))
    stop("reserved character (escape sequence '\\a') detected in expression")
  # Identify replaceable identifiers
  tmp= gregexpr(pattern=rodeoConst$identifierPatterns$core,text=expr)[[1]]
  if (tmp[1] == -1) {
    return(expr) # nothing to substitute
  } else {
    pos= data.frame(stringsAsFactors=FALSE,
      first= tmp,
      final= tmp-1+attr(tmp,which="match.length",exact=TRUE)
    )
    ident= substring(expr, pos$first, pos$final)
    ident= unique(ident)
    # more identifiers than substitutes ?
    bad= ident[!(ident %in% names(sub))]
    if (all && (length(bad) > 0)) {
      stop(paste0("missing substitute(s) for identifier(s) '",
        paste(bad,collapse="', '"),"'"))
    }
    ident= ident[ident %in% names(sub)]
    # Substitute
    # We do this in two steps to avoid the case where (part of) an inserted
    # substitute is later replaced by another substitute
    if (length(ident) > 0) {
      for (i in 1:length(ident)) {
        expr= gsub(pattern=paste0(rodeoConst$identifierPatterns$before,ident[i],
          rodeoConst$identifierPatterns$after),
          replacement=paste0("\\1",specialChar,i,specialChar,"\\2"),
          x=expr)
      }
      for (i in 1:length(ident)) {
        expr= gsub(pattern=paste0(specialChar,i,specialChar),
          replacement=paste0(sub[ident[i]]),
          x=expr)
      }
    }
    return(expr)
  }
}

#substituteIdentifiers("a + b *cd", sub=c(a="a1", b="b2", cd="99", zz="1"), all=TRUE)
#substituteIdentifiers("a + b *cd", sub=c(a="a1", b="b2", cd="99"), all=TRUE)
#substituteIdentifiers("a + b *cd", sub=c(a="a1", b="b2"), all=FALSE)
#substituteIdentifiers("a + b *cd", sub=c(a="a1", b="b2"), all=TRUE)

# Language specific code elements
codeElem= function(lang) {
  if (lang == rodeoConst$lang["r"]) {
    return( list(com="#", cont="", eleOpen="[", eleClose="]",
      vecOpen="c(", vecClose=")", listElem="$", min="min", max="max") )
  } else if (lang == rodeoConst$lang["fortran"]) {
    return( list(com="!", cont="&", eleOpen="(", eleClose=")",
      vecOpen="(/", vecClose="/)", listElem="%", min="min", max="max") )
  } else {
    stop(paste0("target language '",lang,"' not supported; must be one of: '",
      paste(rodeoConst$lang, collapse="', '"),"'"))
  }
}

# Break long Fortran lines
fortran.breakLine= function(text, conti, newline) {
  minlen= 60
  buf=""
  from=1
  k= 0
  for (i in 1:nchar(text)) {
    k= k+1
    if (substr(text,i,i) %in% c("+","-","*","/",",") && (k >= minlen)) {
      if (substr(text,i,min(i+1, nchar(text))) != "**") {
        k= 0
        buf= paste0(buf,substr(text,from,i),conti,newline)
        from=i+1
      }
    }
  }
  if (from <= nchar(text))
    buf= paste0(buf,substr(text,from,nchar(text)))
  return(buf) 
}

# Convert numeric constants into valid Fortran double precision constants
# Notes: This converts both real and integer constants to double precision.
#   It is important to prevent integer divisions) or loss-of-precision problems.
#   See the following test code for supported notations of numeric constants.
##numbers= c("1", "-1", "1e5", "1e-05", "1.", "1.0", "1.0e0", "-1.0e+0", ".1", ".1e0", ".1e+0")
##numbers= paste("prefix99 **",numbers)
##numbers= paste(numbers, " / suffix")
##print(fortran.doubleConst(numbers))
fortran.doubleConst= function(text) {
  # Step 1: Identify numeric constants and enclose within angle brackets
  before= "(^|[^a-zA-Z0-9_])"
  after= "([^a-zA-Z0-9_]|$)"
  pattern= paste0(before,"((?:(?:[-]?[0-9]+[.]?[0-9]*)|(?:[-]?[.][0-9]+))(?:e[-+]?[0-9]+)?)",after)
  replace= "\\1<\\2>\\3"
  text= gsub(pattern=pattern, replacement=replace, x=text)
  # Step 2: Replace exponent symbol "e" by "d"
  pattern= "([<][^>]+)([e])([^<]+[>])"
  replace= "\\1d\\3"
  text= gsub(pattern=pattern, replacement=replace, x=text)
  # Step 3: Append "d0" to constants not written in exponent form
  pattern= "([<])([^d<]+)([>])"
  replace= "\\1\\2d0\\3"
  text= gsub(pattern=pattern, replacement=replace, x=text)
  # Step 4: Strip angle brackets
  pattern= "[<>]"
  replace= ""
  text= gsub(pattern=pattern, replacement=replace, x=text)
  return(text)
}

# Convert power operator ^ into **
fortran.powerOperator= function(text) {
  return(gsub(pattern="^", replacement="**", x=text, fixed=TRUE))
}

