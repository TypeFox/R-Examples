### R code from vignette source 'ex-tcltk-pack.Rnw'

###################################################
### code chunk number 1: ex-tcltk-pack.Rnw:1-4
###################################################
library(tcltk)
## pack examples
## how to pack in buttons


###################################################
### code chunk number 2: ex-tcltk-pack.Rnw:7-11
###################################################
window <- tktoplevel()
tkwm.title(window, "Examples using pack as a geometry manager")
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 3: ex-tcltk-pack.Rnw:31-36
###################################################
frame_1 <- ttklabelframe(frame, text="plain vanilla")
tkpack(frame_1, expand = TRUE, fill = "x")
l <- function(f) 
  list(ttkbutton(f, text="cancel"), ttkbutton(f, text="ok"))
sapply(l(frame_1), tkpack, side = "left")


###################################################
### code chunk number 4: moveRight
###################################################
frame_2 <- ttklabelframe(frame, text = "push to right")
tkpack(frame_2, expand = TRUE, fill = "x")
tkpack(ttklabel(frame_2, text = " "), 
       expand = TRUE, fill = "x", side = "left")
sapply(l(frame_2), tkpack, side = "left")


###################################################
### code chunk number 5: appleSays
###################################################
frame_3 <- ttklabelframe(frame, text="push right with space")
tkpack(frame_3, expand = TRUE, fill = "x")
tkpack(ttklabel(frame_3, text = " "), expand=TRUE, fill="x", 
       side = "left")
sapply(l(frame_3), tkpack, side = "left", padx = 6)


