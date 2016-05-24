

###################################################
### code chunk number 138: Controls.Rnw:1554-1565
###################################################
window <- gwindow("Popup example")
button <- gbutton("click me or right click me", cont = window, 
                  handler = function(h, ...) {
                    cat("You clicked me\n")
                  })
f <- function(h,...) cat("you right clicked on", h$action, "\n")
menu_bar_list <- 
  list(one = gaction("one", action = "one", handler = f),
       two = gaction("two", action = "two", handler = f)
       )
add3rdMousePopupmenu(button, menu_bar_list)
