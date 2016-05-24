
test.gimage <- function() {
  require(RUnit)
  
  w <- gwindow()
  g <- ggroup(cont = w, horizontal=FALSE)

  ## stock icon
  im <- gimage("help", dirname="stock", cont = g)

  checkEquals(svalue(im), system.file("images/help.gif",package="gWidgets"))

  ## svalue<-
  svalue(im) <- "open"

  ## size
  im <- gimage("help", dirname="stock", size="menu", cont = g)
  im <- gimage("help", dirname="stock", size="small_toolbar", cont = g)
  im <- gimage("help", dirname="stock", size="large_toolbar", cont = g)
  im <- gimage("help", dirname="stock", size="button", cont = g)
  im <- gimage("help", dirname="stock", size="dialog", cont = g)
  
  ## filenames
  im <- gimage("test.png", dirname=".", cont = g)

  ## svalue<-
  svalue(im) <- "/Users/verzani/export/Statistics/R/pmg/pmg3/projects/testa.png"
}

test.icons <- function() {
  
  icon <- system.file("images/help.gif",package="gWidgets")

  ## addStockIcons
  addStockIcons("testing",icon)
  im <- gimage("testing", dirname="stock", cont = T)
  checkEquals(svalue(im), icon)

  ## getStockIcons
  lst <- getStockIcons()
  checkEquals(lst[['testing',exact=TRUE]], icon)
  
}
