## RUnit test
## prefix for function is test.

test.gcheckboxgroup <- function() {
  w <- gwindow()
  g <- ggroup(cont = w, horizontal = FALSE)

  items <- letters[1:5]
  
  cbg <- gcheckboxgroup(items, checked=c(TRUE,TRUE, FALSE,FALSE,FALSE), cont = g)

  ## svalue
  checkEquals(svalue(cbg), items[1:2])
  checkEqualsNumeric(svalue(cbg, index=TRUE), 1:2)

  ## svalue<-

  ## by logical
  svalue(cbg) <- c(TRUE,rep(FALSE,4))
  checkEquals(svalue(cbg), items[1])

  ## by name
  svalue(cbg) <- items[1:2]
  checkEquals(svalue(cbg), items[1:2])

  ## by index
  svalue(cbg, index=TRUE) <- 1:3
  checkEquals(svalue(cbg), items[1:3])

  ## [
  checkEquals(cbg[1:2], items[1:2])
  checkEquals(cbg[], items[])

  ## [<-
  items <- toupper(items)
  cbg[] <- items
  checkEquals(cbg[], items[])

  ## clean up
  dispose(w)
}
