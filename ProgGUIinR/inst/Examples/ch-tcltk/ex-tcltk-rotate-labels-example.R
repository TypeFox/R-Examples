
###################################################
### code chunk number 64: rotateLabels
###################################################
window <- tktoplevel()
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both")
#
x <- strsplit("Lorem ipsum dolor sit amet ...", "\\s")[[1]]
labels <- lapply(x, function(i) ttklabel(frame, text = i))
sapply(labels, function(i) tkpack(i, side = "left"))
#
rotateLabel <- function() {
  children <- as.character(tkpack("slaves", frame))
  tkpack.forget(children[1])
  tkpack(children[1], after = children[length(children)], 
         side = "left")
}


###################################################
### code chunk number 65: BasicContainers.Rnw:341-342 (eval = FALSE)
###################################################
for(i in 1:20) {rotateLabel(); Sys.sleep(1)}
