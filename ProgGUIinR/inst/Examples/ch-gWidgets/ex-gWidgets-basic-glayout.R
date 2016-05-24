
###################################################
### code chunk number 40: glayoutExample
###################################################
window <- gwindow("glayout example", visible = FALSE)
lyt <- glayout(cont = window, spacing = 5)
right <- c(1,0); left <- c(-1,0)
lyt[1,1, anchor = right] <- "name"
lyt[1,2, anchor = left ] <- gedit("George Washington", 
           cont = lyt)
#
lyt[2,1, anchor = right] <- "rank"
lyt[2,2, anchor = left ] <- gedit("General", cont = lyt)
#
lyt[3,1, anchor = right] <- "serial number"
lyt[3,2, anchor = left ] <- gedit("1", cont = lyt)
visible(window) <- TRUE


###################################################
### code chunk number 41: main_table_prop
###################################################
sapply(lyt[,2], svalue)
