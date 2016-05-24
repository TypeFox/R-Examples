## From: Christoph Buser <buser@stat.....ethz.ch>
## To: maechler@....
## Subject: Duplicated
## Date: Tue, 25 Sep 2007 14:29:46 +0200

### Changes and more arguments: entirely by MM
Duplicated <- function(v, incomparables = FALSE, fromLast = FALSE,
                       nomatch = NA_integer_)
{
  ## Purpose: A counting-generalization of duplicated()
  ## ----------------------------------------------------------------------
  ## Arguments: a numeric vector
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler & Christoph Buser, Date: Sep 2007

  uv <- unique(nv <- na.omit(v))
  ## easier (but less general?): uv <- unique(nv <- v[!is.na(v)])
  fv <- factor(nv, levels = uv)
  dup <- duplicated(as.integer(fv),
		    incomparables = incomparables, fromLast = fromLast)
  match(v, nv[dup], incomparables = incomparables, nomatch = nomatch)
}
