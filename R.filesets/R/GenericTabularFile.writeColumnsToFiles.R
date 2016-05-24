###########################################################################/**
# @set "class=GenericTabularFile"
# @RdocMethod writeColumnsToFiles
#
# @title "Read each column from a data file and exports it to a separate file"
#
# \description{
#  @get "title".
#  Since each column is processed independently of the others, this method
#  is memory efficient and can handle very large data files.
# }
#
# @synopsis
#
# \arguments{
#   \item{destPath}{The output directory where to write the files.}
#   \item{filenameFmt}{An @see "base::sprintf" format string used to generate
#    filenames given the fullnames (column names plus tags).}
#   \item{tags}{An optional @character @vector of tags added to the fullnames.}
#   \item{columnName}{...}
#   \item{header}{An optional file header.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns (invisibly) a @character @vector of all output files.
# }
#
# \details{
#  Each file generated is written atomically by first writing to a temporary
#  file which is then renamed if successfully written.  This minimizes the
#  risk for creating incomplete files, which otherwise may occur if for
#  instance an interrupt occured.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("writeColumnsToFiles", "GenericTabularFile", function(this, destPath, filenameFmt="%s.txt", tags=NULL, columnName=NULL, header=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  writeHeaderComments0 <- function(con, hdr, commentPrefix="# ", ...) {
    hdr <- c(list(nbrOfHeaderRows=length(hdr)+1, hdr));
    hdrStr <- unlist(hdr);
    hdrStr <- paste(names(hdrStr), hdrStr, sep="\t");
    hdrStr <- paste(commentPrefix, hdrStr, sep="");
    writeLines(con=con, hdrStr);
  }

  escapeFilename <- function(filename, ...) {
    filename <- gsub(":", "%3A", filename);
    filename <- gsub(";", "%3B", filename);
    filename;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'destPath':
  destPath <- Arguments$getWritablePath(destPath);

  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- unlist(strsplit(tags, split=","));
    tags <- trim(tags);
    tags <- tags[nchar(tags, type="chars") > 0L];
  }

  # Argument 'filenameFmt':
  filenameFmt <- Arguments$getCharacter(filenameFmt);

  # Argument 'columnName':
  if (!is.null(columnName))
    columnName <- Arguments$getCharacter(columnName);

  # Argument 'header':
  if (is.null(header)) {
    header <- list(
      sourceFile=getFilename(this)
    );
  } else {
    header <- as.list(header);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  hdrColumnName <- columnName;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify column names
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  columnNames <- getColumnNames(this);
  verbose && printf(verbose, "Column names [%d]:\n", length(columnNames));
  verbose && print(verbose, columnNames);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract and export each column
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathnames <- c();

  colClasses <- "character";
  for (cc in seq_along(columnNames)) {
    columnName <- columnNames[cc];
    verbose && enter(verbose, sprintf("Column #%d ('%s') of %d",
                                       cc, columnName, length(columnNames)));

    fullname <- paste(c(columnName, tags), collapse=",");
    filename <- sprintf(filenameFmt, fullname);
    filename <- escapeFilename(filename);
    pathname <- file.path(destPath, filename);

    # Check if file already exists
    if (!isFile(pathname)) {
      names(colClasses) <- sprintf("^%s$", columnName);
      values <- readDataFrame(this, colClasses=colClasses);
      values <- trim(values[[1]]);
      df <- data.frame(dummy=values, stringsAsFactors=FALSE);
      if (is.null(hdrColumnName)) {
        colnames(df) <- columnName;
      } else {
        colnames(df) <- hdrColumnName;
      }
      verbose && str(verbose, df);

      # Write atomically, by first writing to a temporary file
      pathnameT <- sprintf("%s.tmp", pathname);
      pathnameT <- Arguments$getWritablePathname(pathnameT, mustNotExist=TRUE);

      con <- file(pathnameT, open="w");
      header$createdOn <- format(Sys.time(), "%Y-%m-%d %H:%M:%S");
      header$column <- cc;
      header$columnName <- columnName;
      header$nbrOfDataRows <- nrow(df);
      writeHeaderComments0(con=con, header);
      write.table(file=con, df, quote=FALSE, sep="\t", row.names=FALSE);
      close(con);

      # Rename temporary file
      verbose && enter(verbose, "Renaming temporary file to destination name");
      res <- file.rename(pathnameT, pathname);
      if (!res) {
        throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
      }
      verbose && exit(verbose);
    } else {
      verbose && cat(verbose, "Column already extracted");
    }

    pathnames <- c(pathnames, pathname);

    verbose && exit(verbose);
  } # for (cc ...)

  invisible(pathnames);
})


############################################################################
# HISTORY:
# 2012-12-03
# o Generalized writeColumnsToFiles() to GenericTabularFile.  Used to
#   be only for TabularTextFile.
# 2011-02-18
# o ROBUSTNESS: Now writeColumnsToFiles() for TabularTextFile writes
#   files atomically, which should minimize the risk for generating
#   incomplete files.
# 2009-04-17
# o Added Rdoc comments.
# 2008-05-21
# o BUG FIX: Argument 'verbose' was never passed to Arguments$getVerbose().
# 2008-05-05
# o Now some non-valid filename characters are escaped.
# o Added internal escapeFilename().
# 2008-05-01
# o Created.
############################################################################
