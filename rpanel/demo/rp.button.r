rp.messagebox("Default pos")

p1 <- rp.control()
rp.button(p1, function(p1){ rp.messagebox("Pressed!"); p1 }, "Press me")

rp.messagebox("Pos is Left")

p2 <- rp.control()
rp.button(p2, function(p2){ rp.messagebox("Pressed!"); p2 }, "Press me", pos="left")

rp.messagebox("Pos is Right")

p3 <- rp.control()
rp.button(p3, function(p3){ rp.messagebox("Pressed!"); p3 }, "Press me", pos="right")

rp.messagebox("Pos is x, y, w, h")

p4 <- rp.control()
rp.button(p4, function(p4){ rp.messagebox("Pressed!"); p4 }, "Press me", pos=list(100,100,200,50))

rp.messagebox("Pos is grid")

p5 <- rp.control(a=10)
rp.button(p5, function(p5){ rp.messagebox("Pressed!"); p5 }, "Messagebox", pos=list(row=0, column=0, sticky="news"), name="b1")
rp.button(p5, function(p5){ p5$a <- 10; p5 }, "Set to 10", pos=list(row=1, column=1, sticky="s"), background="Red", foreground="White", font="Arial", name="b2")
rp.button(p5, function(p5){ p5$a <- p5$a+1; print(p5$a); p5 }, "Add 1", pos=list(row=2, column=2), name="b3")
rp.button(p5, function(p5){ print(p5$a); p5 }, "Print value", pos=list(row=3, column=0), , name="b4")
rp.button(p5, function(p5){ rp.widget.dispose(p5, "b1"); rp.widget.dispose(p5, "b2"); rp.widget.dispose(p5, "b3"); rp.widget.dispose(p5, "b4"); p5 }, "Remove buttons", pos=list(row=4, column=0))
rp.button(p5, function(p5){ rp.control.dispose(p5); p5 }, "Close window", pos=list(row=5, column=0))

