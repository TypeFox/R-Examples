
###################################################
### code chunk number 37: Containers.Rnw:47-52
###################################################
window <- gtkWindow(show=FALSE)          # use default type
window$setTitle("Window title")          # set window title
window['title']                          # or  use getTitle
window$setDefaultSize(250,300)           # 250 wide, 300 high
window$show()                            # show window
