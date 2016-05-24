
###################################################
### code chunk number 45: Overview.Rnw:1192-1201
###################################################
window <- tktoplevel()
button <- 
  ttkbutton(window, text = "Click me for the x,y position")
tkpack(button)
tkbind(button, "<ButtonPress-1>", function(W, x, y, X, Y) {
  print(W)                              # an ID
  print(c(x, X))                        # character class
  print(c(y, Y))
  })
