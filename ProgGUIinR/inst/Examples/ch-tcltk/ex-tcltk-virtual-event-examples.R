
###################################################
### code chunk number 42: bind_examples
###################################################
window <- tktoplevel()
label <- ttklabel(window, text = "Click Ok for a message")
button1 <- ttkbutton(window, text = "Cancel", 
                command = function() tkdestroy(window))
button2 <- ttkbutton(window, text = "Ok", command=function() {
  print("initiate an action")
})
sapply(list(label, button1, button2), tkpack)
##
tkbind(window, "C", function() tcl(button1, "invoke"))
tkconfigure(button1, underline = 0)
##
tkbind(window, "O", function() tcl(button1, "invoke"))
tkconfigure(button2, underline = 0)
tkfocus(button2)
##
tkbind("TButton", "<Return>", function(W) {
  tcl(W, "invoke")
})


###################################################
### code chunk number 43: Overview.Rnw:1112-1114
###################################################
 tkevent.add("<<Paste>>", "<Control-y>")
 tkevent.add("<<Save>>", "<Control-x><Control-s>")

