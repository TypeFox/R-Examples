setConstructorS3("RichDataFrame", function(..., .name=NULL, .virtuals=list()) {
  data <- data.frame(...);
  this <- extend(data, "RichDataFrame");
  this <- setName(this, .name);
  this <- setVirtualColumnFunctions(this, .virtuals);
  this;
})


setMethodS3("getGenericSummary", "RichDataFrame", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1L]);
  name <- getName(this);
  if (is.null(name)) name <- "<none>";
  s <- c(s, sprintf("Name: %s", as.character(name)));
  s <- c(s, sprintf("Number of rows: %d", nrow(this)));

  data <- this[integer(0L),,drop=FALSE];
  data <- as.data.frame(data);
  names <- getColumnNames(this);
  vNames <- getVirtualColumnNames(this);
  names <- sapply(names, FUN=function(field) {
    values <- data[[field]];
    mode <- mode(values);
    if (is.element(field, vNames)) {
      field <- sprintf("%s*", field);
    }
    sprintf("%s [%s]", field, mode);
  });
  if (length(names) == 0L) names <- "<none>";

  s <- c(s, sprintf("Number of columns: %d", length(names)));
  names <- paste(names, collapse=", ");
  s <- c(s, sprintf("Columns: %s", names));

  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  GenericSummary(s);
}, private=TRUE)


setMethodS3("print", "RichDataFrame", function(x, ...) {
  print(as.data.frame(x, ...));
})

setMethodS3("setColumnNamesTranslator", "RichDataFrame", function(this, fcn, ...) {
  # Argument 'fcn':
  if (!is.null(fcn)) {
    stopifnot(is.function(fcn));
  }

  attr(this, ".columnNamesTranslator") <- fcn;
  invisible(this);
})

setMethodS3("setColumnNamesMap", "RichDataFrame", function(this, ...) {
  map <- c(...);

  mapInv <- names(map);
  names(mapInv) <- map;

  fcn <- function(names, invert=FALSE, ...) {
    n <- length(names);

    # Nothing to do?
    if (n == 0) {
      return(names);
    }

    # Invert?
    if (invert) {
      map <- mapInv;
    }

    idxs <- match(names, names(map));
    toUpdateIdxs <- which(!is.na(idxs));

    # Nothing to do?
    if (length(toUpdateIdxs) == 0) {
      return(names);
    }

    names[toUpdateIdxs] <- map[names[toUpdateIdxs]];

    # Sanity check
    stopifnot(length(names) == n);

    names;
  } # fcn()

  setColumnNamesTranslator(this, fcn, ...);
}) # setColumnNamesMap()


setMethodS3("getColumnNamesTranslator", "RichDataFrame", function(this, ...) {
  attr(this, ".columnNamesTranslator");
})


setMethodS3("translateColumnNames", "RichDataFrame", function(this, names, ...) {
  fcn <- getColumnNamesTranslator(this);
  if (!is.null(fcn)) {
    names <- fcn(names, ...);
  }
  names;
})

setMethodS3("getColumnNames", "RichDataFrame", function(this, virtual=TRUE, translate=TRUE, vTag=NULL, ...) {
  names <- names(unclass(this));
  if (translate) {
    names <- translateColumnNames(this, names);
  }

  if (virtual) {
    vNames <- getVirtualColumnNames(this, translate=translate);
    if (!is.null(vTag)) {
      vNames <- sprintf("%s%s", vNames, vTag);
    }
    names <- c(names, vNames);
  }

  names;
}) # getColumnNames()


setMethodS3("names", "RichDataFrame", function(x) {
  getColumnNames(x);
}, appendVarArgs=FALSE)


setMethodS3("hasColumns", "RichDataFrame", function(this, names, ...) {
  # Argument 'names':
  names <- Arguments$getCharacters(names);
  is.element(names, getColumnNames(this, ...));
}) # hasColumns()

setMethodS3("hasColumn", "RichDataFrame", function(this, name, ...) {
  # Argument 'name':
  name <- Arguments$getCharacter(name);
  hasColumns(this, name);
}) # hasColumn()



setMethodS3("getFullName", "RichDataFrame", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this, collapse=NULL);
  paste(c(name, tags), collapse=",");
})


setMethodS3("getName", "RichDataFrame", function(this, ...) {
  attr(this, "name", exact=TRUE);
})


setMethodS3("getTags", "RichDataFrame", function(this, collapse=",", ...) {
  tags <- attr(this, "tags", exact=TRUE);
  if (length(tags) == 0) {
    return(tags);
  }
  tags <- unlist(strsplit(tags, split=","));
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  }
  tags;
})


setMethodS3("setName", "RichDataFrame", function(this, name, ...) {
  if (!is.null(name)) {
    name <- Arguments$getCharacter(name);
  }
  attr(this, "name") <- name;
  invisible(this);
})


setMethodS3("setTags", "RichDataFrame", function(this, tags, ...) {
  # Argument 'tags':
  tags <- Arguments$getCharacters(tags);
  tags <- unlist(strsplit(tags, split=","));
  tags <- tags[nchar(tags) > 0];
  attr(this, "tags") <- tags;
  invisible(this);
})


setMethodS3("setAttributes", "RichDataFrame", function(this, attrs=NULL, ...) {
  # Argument 'attrs':
  if (is.null(attrs)) {
    return(invisible(this));
  }
  stopifnot(is.list(attrs));

  exclNames <- c("class", "names", "row.names");
  names <- names(attrs);
  names <- setdiff(names, exclNames);
  attrs <- attrs[names];
  for (key in names(attrs)) {
    attr(this, key) <- attrs[[key]];
  }
  invisible(this);
})


setMethodS3("newInstance", "RichDataFrame", function(this, nrow, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'nrow':
  nrow <- Arguments$getInteger(nrow, range=c(0,Inf));

  # Allocate empty data frame populated with missing values
  res <- this;
  res[1L,] <- NA;
  subset <- rep(1L, times=nrow);
  res <- res[subset,,drop=FALSE];
  rownames(res) <- NULL;

  # Preserve all attributes
  setAttributes(res, attributes(this));

  # Sanity check
  stopifnot(nrow(res) == nrow);

  res;
}, protected=TRUE) # newInstance()


setMethodS3("subset", "RichDataFrame", function(x, subset, select, drop=FALSE, envir=parent.frame(), ...) {
  # To please R CMD check
  this <- x;

  # Data table to use for indentification of rows and columns
  data <- as.data.frame(this);

  # Identify rows to select
  if (missing(subset)) {
    r <- TRUE;
  } else {
    e <- substitute(subset);
    r <- eval(e, envir=data, enclos=envir);
    if (!is.logical(r)) {
      throw("Argument 'subset' must evaluate to logical: ", mode(r));
    }
    r <- r & !is.na(r);
  }

  # Identify columns to select
  if (missing(select)) {
    vars <- TRUE;
  } else {
    nl <- as.list(seq_along(data));
    names(nl) <- names(data);
    vars <- eval(substitute(select), envir=nl, enclos=envir);
  }

  # Not needed anymore
  data <- NULL;

  res <- this[r,, drop=drop];

  res;
})


setMethodS3("as.data.frame", "RichDataFrame", function(x, ..., virtual=TRUE, translate=TRUE) {
  # To please R CMD check
  this <- x;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return "default" columns
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- this;
  attr(data, ".virtuals") <- NULL;
  class(data) <- "data.frame";

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Append virtual columns
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (virtual) {
    vNames <- getVirtualColumnNames(this);
    if (length(vNames) > 0) {
      dataT <- data;

      # Translate column names?
      names0 <- names(data);
      names <- translateColumnNames(this, names0);
      if (!all(names == names0)) {
        names(dataT) <- names;
      }

      for (name in vNames) {
        data[[name]] <- getVirtualColumn(this, name, data=dataT);
      } # for (name ...)

      # Not needed anymore
      dataT <- NULL;
    }
  }

  if (translate) {
    names <- names(data);
    names <- translateColumnNames(this, names);
    names(data) <- names;
  }

  data;
})

setMethodS3("dim", "RichDataFrame", function(x) {
  names <- names(x);
  ncol <- length(names);
  nrow <- length(unclass(x)[[1]]);
  c(nrow, ncol);
}, appendVarArgs=FALSE)

setMethodS3("length", "RichDataFrame", function(x) {
  names <- names(x);
  length(names);
}, appendVarArgs=FALSE)

setMethodS3("as.list", "RichDataFrame", function(x, ...) {
  x <- as.data.frame(x);
  as.list(x, ...);
})


setMethodS3("[", "RichDataFrame", function(x, i, j, drop=NULL) {
  this <- x;

  virtuals <- getVirtualColumnFunctions(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Columns
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  colnames <- getColumnNames(this);
  ncol <- length(colnames);
  cols <- seq_len(ncol);
  names(cols) <- colnames;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Rows
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nrow <- nrow(this);
  rows <- seq_len(nrow);
  names(rows) <- rownames(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Subset rows and columns
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Dummy matrix to catch out of bound errors
  data <- matrix(NA, nrow=nrow, ncol=ncol);
  colnames(data) <- names(cols);
  rownames(data) <- names(rows);

  if (missing(i) && missing(j)) {
  } else if (missing(i)) {
    data <- data[,j,drop=FALSE];
    cols <- cols[j];
  } else if (missing(j)) {
    data <- data[i,,drop=FALSE];
    rows <- rows[i];
  } else {
    data <- data[i,j,drop=FALSE];
    rows <- rows[i];
    cols <- cols[j];
  }

  # Not needed anymore
  data <- NULL;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Drop?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(drop)) {
    drop <- FALSE;
    if (missing(i)) {
      drop <- TRUE;
    } else {
      drop <- (length(cols) == 1);
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Virtual columns
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  vColnames <- getVirtualColumnNames(this);
  nvcol <- length(vColnames);

  # Are there any virtual columns?
  if (nvcol > 0) {
    isVirtual <- rep(FALSE, times=ncol);
    isVirtual[is.element(colnames, vColnames)] <- TRUE;
    names(isVirtual) <- colnames;

    # Subset virtual columns
    vColnamesS <- intersect(vColnames, colnames[cols]);
    virtuals <- virtuals[vColnamesS];

    # Subset non-virtual columns
    cols <- cols[!isVirtual[cols]];
  } # if (nvcol > 0)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Subset
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  i <- rows;
  j <- cols;
##  res <- NextMethod("[");

  class <- class(this);
  res <- this;
  class(res) <- "data.frame";
  res <- res[i,j,drop=drop];


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Preserve attributes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Keep class attributes, unless dimensions have been dropped
  if (!is.null(dim(res))) {
    class(res) <- class;
    res <- setAttributes(res, attributes(this));
    res <- setVirtualColumnFunctions(res, virtuals);
  }

  res;
})


setMethodS3("[[", "RichDataFrame", function(x, name) {
  this <- x;

  colnames <- getColumnNames(this);
  ncol <- length(colnames);

  # Argument 'name':
  if (length(name) == 0) {
    throw("attempt to select less than one element");
  }
  if (length(name) > 1) {
    throw("attempt to select more than one element");
  }
  if (is.numeric(name)) {
    if (name < 1 || name > ncol) {
      throw(sprintf("Column index out of range [%d,%d]: %d", min(1,ncol), ncol, name));
    }
    name <- colnames[name];
  } else {
    if (!is.element(name, colnames)) {
      return(NULL);
    }
  }

  # Is it a virtual column?
  vColnames <- getVirtualColumnNames(this);
  isVirtual <- is.element(name, vColnames);
  if (isVirtual) {
    res <- getVirtualColumn(this, name);
  } else {
    col <- match(name, colnames);
    res <- unclass(this)[[col]];
  }

  res;
})

setMethodS3("[[<-", "RichDataFrame", function(x, name, value) {
  this <- x;

  name <- translateColumnNames(this, name, invert=TRUE);

  # Setting a virtual column?
  if (is.function(value)) {
    if (is.element(name, getColumnNames(this, virtual=FALSE, translate=FALSE))) {
      this[[name]] <- NULL;
    }
    this <- setVirtualColumn(this, name, value);
  } else {
    # Drop existing virtual column?
    if (hasVirtualColumn(this, name, translate=FALSE)) {
      this <- dropVirtualColumn(this, name);
    }

    # Assign & preserve attributes
    attrs <- attributes(this);
    tryCatch({
      this <- NextMethod("[[<-");
    }, error=function(ex) {
      msg <- sprintf("Failed assign column '%s', because: %s", name, ex$message);
      throw(msg);
    })
    this <- setAttributes(this, attrs);
  }

  invisible(this);
})


setMethodS3("$", "RichDataFrame", function(x, name) {
  x[[name]];
})

setMethodS3("$<-", "RichDataFrame", function(x, name, value) {
  x[[name]] <- value;
  invisible(x);
})



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# VIRTUAL COLUMNS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("getVirtualColumnFunctions", "RichDataFrame", function(this, translate=TRUE, ...) {
  res <- attr(this, ".virtuals", exact=TRUE);
  if (is.null(res)) {
    res <- list();
  }
  # Sanity check
  stopifnot(is.list(res));

  if (translate) {
    names <- names(res);
    names <- translateColumnNames(this, names);
    names(res) <- names;
  }

  res;
}, private=TRUE) # getVirtualColumnFunctions()


setMethodS3("getVirtualColumnFunction", "RichDataFrame", function(this, name, ...) {
  # Argument 'name':
  name <- Arguments$getCharacter(name);

  virtuals <- getVirtualColumnFunctions(this, ...);
  if (!is.element(name, names(virtuals))) {
    throw("No such virtual column: ", name);
  }

  fcn <- virtuals[[name]];

  # Sanity check
  stopifnot(is.function(fcn));

  fcn;
})



setMethodS3("setVirtualColumnFunctions", "RichDataFrame", function(this, virtuals, ...) {
  # Argument 'virtuals':
  if (is.null(virtuals)) {
    return(invisible(this));
  }

  stopifnot(is.list(virtuals));

  if (length(virtuals) > 0) {
    # Must be named
    stopifnot(!is.null(names(virtuals)));

    # Cannot have duplicated names
    stopifnot(!any(duplicated(names(virtuals))));

    # Assert non-name clashes between virtual and default columns.
    stopifnot(all(!is.element(names(virtuals), getColumnNames(this, virtual=FALSE))));
  }

  # Validate content
  for (key in names(virtuals)) {
    fcn <- virtuals[[key]];
    if (is.function(fcn)) {
    } else if (is.character(fcn)) {
      fcn <- Arguments$getCharacter(fcn);
    } else {
      throw(sprintf("Unknown type of virtual column '%s': %s", key, mode(fcn)));
    }
  } # for (kk ...)

  attr(this, ".virtuals") <- virtuals;

  invisible(this);
}, private=TRUE) # setVirtualColumnFunctions()


setMethodS3("getVirtualColumnNames", "RichDataFrame", function(this, ...) {
  virtuals <- getVirtualColumnFunctions(this, ...);
  names(virtuals);
}) # getVirtualColumnNames()


setMethodS3("hasVirtualColumns", "RichDataFrame", function(this, names, ...) {
  # Argument 'names':
  names <- Arguments$getCharacters(names);

  virtuals <- getVirtualColumnFunctions(this, ...);
  is.element(names, names(virtuals));
}) # hasVirtualColumnNames()


setMethodS3("hasVirtualColumn", "RichDataFrame", function(this, name, ...) {
  # Argument 'name':
  name <- Arguments$getCharacter(name);

  hasVirtualColumns(this, name, ...);
}) # hasVirtualColumnName()



setMethodS3("getVirtualColumn", "RichDataFrame", function(this, name, data=NULL, ...) {
  fcn <- getVirtualColumnFunction(this, name, ...);

  if (is.null(data)) {
    data <- as.data.frame(this, virtual=FALSE);
  }

  tryCatch({
    value <- fcn(data, ...);
  }, error=function(ex) {
    msg <- sprintf("Failed to calculate virtual column '%s', because: %s", name, ex$message);
    throw(msg);
  })

  # Sanity check
  if (length(value) != nrow(data)) {
    throw(sprintf("The function for virtual column '%s' did not return the expected number of values: %d != %d", name, length(value), nrow(data)));
  }

  value;
}, protected=TRUE)


setMethodS3("dropVirtualColumn", "RichDataFrame", function(this, name, ...) {
  # Argument 'name':
  if (!hasVirtualColumn(this, name, ...)) {
    throw("No such virtual column: ", name);
  }

  virtuals <- getVirtualColumnFunctions(this);
  virtuals[[name]] <- NULL;
  this <- setVirtualColumnFunctions(this, virtuals);

  invisible(this);
})


setMethodS3("setVirtualColumn", "RichDataFrame", function(this, name, fcn, ...) {
  # Argument 'name':
  name <- Arguments$getCharacter(name);

  # Argument 'fcn':
  if (is.character(fcn)) {
    if (!exists(fcn, mode="function")) {
      throw(sprintf("Cannot set virtual column '%s'. No such function: %s", name, fcn));
    }
    fcn <- get(fcn, mode="function");
  }

  if (!is.function(fcn)) {
    throw(sprintf("Cannot set virtual column '%s'. Argument 'fcn' is not specifying a function: ", mode(fcn)));
  }

  virtuals <- getVirtualColumnFunctions(this);
  virtuals[[name]] <- fcn;
  this <- setVirtualColumnFunctions(this, virtuals);

  invisible(this);
}, protected=TRUE)


setMethodS3("rbind", "RichDataFrame", function(..., deparse.level=1) {
  # Argument '...':
  args <- list(...);
  this <- args[[1]];
  class <- class(this);
  for (kk in seq_along(args)) {
    other <- args[[kk]];
    other <- Arguments$getInstanceOf(other, class[1]);
  } # for (kk ...)

  # Nothing to do?
  if (length(args) == 1) {
    return(this);
  }

  # Append
  class(this) <- "data.frame";
  class(other) <- "data.frame";
  res <- rbind(this, other);
  class(res) <- class;

  # Preserve attributes
  res <- setAttributes(res, attributes(this));

  res;
}) # rbind()


############################################################################
# HISTORY:
# 2014-01-16
# o BUG FIX: getGenericSummary() for RichDataFrame would throw "Error in
#   data[i, , drop = FALSE] : subscript out of bounds" if it had zero rows.
# 2012-03-23
# o Now argument '.virtuals' of RichDataFrame defaults to list().
# o BUG FIX: "["() for RichDataFrame would loose the class attribute,
#   unless argument 'drop' was FALSE.
# 2012-03-14
# o SPEEDUP: Now "["() no longer uses as.data.frame().
# o SPEEDUP: Now getColumnNames() no longer uses as.data.frame().
# o BUG FIX: "["(i,j) would insert NA rows if any row indices 'i' where
#   duplicated.  Same for columns.
# o BUG FIX: Had to adopt/paste from subset() for data.frame, because
#   there are expressions that are evaluated in the "parent.frame()".
# 2012-03-13
# o Created RichDataFrame from RawGenomicSignals.
############################################################################
