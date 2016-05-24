setMethodS3("getIdentifier", "AffymetrixFileSet", function(this, ...) {
  path <- getPath(this);
  res <- NULL;
  for (kk in 1:3) {
    pathname <- file.path(path, "IDENTIFIER");
    if (isFile(pathname)) {
      res <- readLines(pathname);
      # Remove comments
      res <- trim(gsub("#.*", "", trim(res)));
      # Remove empty lines
      res <- res[nzchar(res)];
      break;
    }
    path <- dirname(path);
  }

  if (!is.null(res)) {
    res <- getChecksum(list(res));
  }

  res;
}, private=TRUE)


############################################################################
# HISTORY:
# 2008-05-10
# o Removed getDescription().
# 2006-11-20
# o Made getDescription() protected.
# 2006-10-30
# o Added getDescription() which search and parse all DESCRIPTION files in
#   the data-set directory tree.
# o Added getIdentifier() which returns a 32-character long hexadecimal
#   hashcode for the "Identifier" string returned by getDescription().
#   If no such string exists, NULL is returned.  This will allow users
#   to specify their own identifiers.
############################################################################
