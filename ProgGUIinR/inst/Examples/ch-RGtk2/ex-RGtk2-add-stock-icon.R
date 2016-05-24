### R code from vignette source 'ex-RGtk2-add-stock-icon.Rnw'

###################################################
### code chunk number 1: ex-RGtk2-add-stock-icon.Rnw:5-10
###################################################
## This example shows, without much explanation the steps to add images
## to the list of stock icons. To generate some sample icons, we use
## those provided by objects in the \pkg{ggplot2} package.

library(RGtk2)


###################################################
### code chunk number 2: ex-RGtk2-add-stock-icon.Rnw:15-34
###################################################
## This example shows how to add new icons to the stock icons
## It uses as an icon source the elements of ggplot2, which have an icon 
## built into them.


## From the API
## An icon factory manages a collection of GtkIconSet; 
## a GtkIconSet manages a set of variants of a particular icon 
## (i.e. a GtkIconSet contains variants for different sizes and 
## widget states). Icons in an icon factory are named by a stock ID, 
## which is a simple string identifying the icon. Each GtkStyle has a 
## list of GtkIconFactory derived from the current theme; those icon 
## factories are consulted first when searching for an icon. If the 
## theme doesn't set a particular icon, GTK+ looks for the icon in a 
## list of default icon factories, maintained by 
## gtk_icon_factory_add_default() and gtk_icon_factory_remove_default(). 
## Applications with icons should add a default icon factory with
## their icons, which will allow themes to override the icons for the 
## application. 


###################################################
### code chunk number 3: makePixbufs
###################################################
require(ggplot2)
iconNames <- c("GeomBar","GeomHistogram")   # 2 of many ggplot functions
icon.size <- 16
pixbufs <- sapply(iconNames, function(name) {
  pixmap <- gdkPixmap(drawable = NULL, width = icon.size, height = icon.size,
                      depth = 24)
  asCairoDevice(pixmap)
  val <- try(get(name))
  grid.newpage()
  try(grid.draw(val$icon()), silent=TRUE)
  dev.off()
  gdkPixbufGetFromDrawable(NULL, pixmap, NULL,0,0,0,0,-1,-1)
})


###################################################
### code chunk number 4: addToStockIcons
###################################################
addToStockIcons <- function(pixbufs, stock.prefix="new") {
  iconfactory <- gtkIconFactory()
  
  items <- lapply(names(pixbufs), function(iconName) {
    ## each icon has its own icon set, which is registered with icon factory
    iconset <- gtkIconSetNewFromPixbuf(pixbufs[[iconName]])
    stockName <- paste(stock.prefix, "-", iconName, sep="")
    iconfactory$add(stockName, iconset)
    
    ## create stock item for icon
    as.GtkStockItem(list(stock_id = stockName, label = iconName))
  })
  ## register our factory of icons
  iconfactory$addDefault()
  ## officially register the stock items
  gtkStockAdd(items)
}


###################################################
### code chunk number 5: ex-RGtk2-add-stock-icon.Rnw:74-77
###################################################
addToStockIcons(pixbufs)
nms <- gtkStockListIds()
unlist(nms[grep("^new", nms)])


