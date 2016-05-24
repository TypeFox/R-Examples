###################################################
### chunk number 1: 
###################################################
someColors <- c("black", "red", "blue", "brown",
                "green", "yellow", "purple",
                paste("grey", seq(10,90,by=10), sep=""))


###################################################
### chunk number 2: 
###################################################
require(grid)
iconDir <- tempdir(); iconSize <- 16;
makeColorIcon <- function(i) {
  filename <- paste(iconDir,"/color-",i,".png",
                    sep="",collapse="")
  png(file=filename,width=iconSize,height=iconSize)
  grid.newpage()
  grid.draw(rectGrob(gp = gpar(fill=i)))
  dev.off()
  return(filename)
}


###################################################
### chunk number 3: 
###################################################
icons <- sapply(someColors, makeColorIcon)
iconNames <- paste("color-",someColors,sep="")
addStockIcons( iconNames, icons)


###################################################
### chunk number 4: 
###################################################
w <- gwindow("icon example")
f <- function(h,...) print(h$action)
tbl <- glayout(cont = w, spacing=0)
for(i in 1:4) {
  for(j in 1:4) {
    ind <- (i - 1) * 4 + j
    tbl[i,j] <- gimage(icons[(i-1)*4 + j], handler=f, 
                       action = iconNames[ind], cont = tbl)
  }
}


