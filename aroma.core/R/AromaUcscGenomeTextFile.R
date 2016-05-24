setConstructorS3("AromaUcscGenomeTextFile", function(...) {
  extend(AromaGenomeTextFile(...), "AromaUcscGenomeTextFile");
})


setMethodS3("findByGenome", "AromaUcscGenomeTextFile", function(static, genome, type, tags=NULL, pattern=sprintf("^%s,UCSC(|,%s),%s[.]txt$", genome, paste(tags, collapse=","), type), ...) {
  NextMethod("findByGenome", genome=genome, tags=NULL, pattern=pattern);
}, static=TRUE)

setMethodS3("readDataFrame", "AromaUcscGenomeTextFile", function(this, colClasses=c("*"=NA), ...) {
  NextMethod("readDataFrame", colClasses=colClasses);
})

############################################################################
# HISTORY:
# 2014-02-03
# o BUG FIX: readDataFrame() for AromaUcscGenomeTextFile explicitly passed
#   arguments '...' to NextMethod(), which would cause them to be
#   duplicated in certain cases.
# 2011-11-17
# o Created.
############################################################################
