setMethodS3("fullname", "default", function(name, tags=NULL, ..., collapse=TRUE) {
  # Create a clean vector of parts
  args <- list(name, tags, ...);
  parts <- unlist(args, use.names=FALSE);
  parts <- paste(parts, collapse=",");
  parts <- strsplit(parts, split=",", fixed=TRUE);
  parts <- unlist(parts, use.names=FALSE);

  parts <- parts[nchar(parts, type="chars") > 0L];

  if (collapse) {
    parts <- paste(parts, collapse=",");
  }

  parts;
})


setMethodS3("name", "default", function(...) {
  # Create a clean vector of parts
  parts <- fullname(..., collapse=FALSE);

  # Extract the name
  name <- parts[1];

  name;
})


setMethodS3("tags", "default", function(..., collapse=FALSE) {
  # Create a clean vector of parts
  parts <- fullname(..., collapse=FALSE);

  # Extract the tags
  tags <- parts[-1];

  if (collapse) {
    tags <- paste(tags, collapse=",");
  }
  tags;
})


setMethodS3("dropTags", "default", function(..., drop=NULL, collapse=FALSE) {
  parts <- fullname(..., collapse=FALSE)

  # Argument 'drop':
  drop <- fullname(drop, collapse=FALSE)

  parts <- c(parts[1], setdiff(parts[-1], drop))

  fullname(parts, collapse=collapse)
})


############################################################################
# HISTORY:
# 2016-01-02
# o BUG FIX: dropTags() would drop name if a tag had the same name.
# 2011-03-09
# o Added dropTags().
# o Now all these functions drops empty tags (and names).
# 2011-03-08
# o Added fullname(), name(), and tags().  Amazing that I never thought
#   of having these very basic functions.  They will make some code
#   and examples much cleaner.  They are also great for illustrating
#   the definition of fullnames, names and tags.
# o Created.
############################################################################
