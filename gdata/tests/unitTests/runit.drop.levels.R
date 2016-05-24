### runit.drop.levels.R
###------------------------------------------------------------------------
### What: Tests for drop.levels
### $Id: runit.drop.levels.R 1784 2014-04-05 02:23:45Z warnes $
### Time-stamp: <2006-08-29 14:21:12 ggorjan>
###------------------------------------------------------------------------

### {{{ --- Test setup ---

if(FALSE) {
  library("RUnit")
  library("gdata")
}

### }}}
### {{{ --- drop.levels ---

test.drop.levels <- function()
{
  f <- factor(c("A", "B", "C", "D"))[1:3]
  fDrop <- factor(c("A", "B", "C"))

  l <- list(f=f, i=1:3, c=c("A", "B", "D"))
  lDrop <- list(f=fDrop, i=1:3, c=c("A", "B", "D"))

  df <- as.data.frame(l)
  dfDrop <- as.data.frame(lDrop)

  checkIdentical(drop.levels(f), fDrop)
  checkIdentical(drop.levels(l), lDrop)
  checkIdentical(drop.levels(df), dfDrop)
}

### }}}
### {{{ Dear Emacs
## Local variables:
## folded-file: t
## End:
### }}}

###------------------------------------------------------------------------
### runit.drop.levels.R ends here
