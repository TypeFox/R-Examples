#example function for "pos" documentation.  Can't be an official example because it calls block and messagebox.

invisible("This function will be called on pressing the buttons.")
invisible("It pops up a message box with the pos parameter inside.")
# This function will be called on pressing the buttons.
# It pops up a message box with the pos parameter inside.
showpos <- function(pos){
  # this function returns the action function: how it 
  # does this is not relevant to the layout demo.
  function(panel,...) {
    rp.messagebox("The position of this button is ",pos,".")
    panel
  }
}

rp.messagebox("This is an interactive demonstration of the use of the 'pos' parameter.")
rp.messagebox("Firstly, this is what happens if you don't specify any 'pos' parameters.  Remember to resize the window and see what happens, and close the window when you're done.")
# Create an rpanel and add some buttons to it: 'default' mode

panel1 <- rp.control(title='Default mode: no pos specified')
rp.button(panel1, action = showpos('NULL'), title = "Button 1")
rp.button(panel1, action = showpos('NULL'), title = "Button 2")
rp.button(panel1, action = showpos('NULL'), title = "Button 3")
rp.block(panel1)

rp.messagebox("Now we see what we can do with the 'pack' options.  Again, note what it does when you resize the window.")

#create an rpanel and put buttons in using "pack" options.
panel2 <- rp.control(title="Pack mode: pos=left/right/top/bottom")
rp.button(panel2, action = showpos("'left'"), title = "Button 1", pos = "left")
rp.button(panel2, action = showpos("'left'"), title = "Button 2", pos = "left")
rp.button(panel2, action = showpos("'right'"), title = "Button 3", pos = "right")
rp.button(panel2, action = showpos("'top'"), title = "Button 4", pos = "top")
rp.button(panel2, action = showpos("'bottom'"), title = "Button 5", pos = "bottom")
rp.button(panel2, action = showpos("'right'"), title = "Button 6", pos = "right")
rp.button(panel2, action = showpos("'bottom'"), title = "Button 7", pos = "bottom")
rp.block(panel2)

rp.messagebox("A slightly more sensible use of this is to have an interactive image or a Tkrplot to the side of some buttons...")

panel3 <- rp.control(title="Graphic and buttons, using 'pack'")
image.file <- file.path(system.file(package = "rpanel"), "images", "gulllmks.gif")
rp.image(panel3, image.file, pos = "left", id = "gulls.image", action = showpos("'left'"))
rp.button(panel3, action = showpos("'top'"), title = "Button 1", pos = "top")
rp.button(panel3, action = showpos("'top'"), title = "Button 2", pos = "top")
rp.button(panel3, action = showpos("'top'"), title = "Button 3", pos = "top")
rp.button(panel3, action = showpos("'top'"), title = "Button 4", pos = "top")
rp.button(panel3, action = showpos("'bottom'"), title = "Button 5", pos = "bottom")
rp.block(panel3)

rp.messagebox("Specifying pixel positions allows a huge amount of control, but takes more time to do and does not resize dynamically")

panel4 <- rp.control(title="'place' demo",size=c(300,150))
rp.button(panel4, action = showpos("c(20,20,40,40"), title = "Button 1", pos = c(20,20,40,40))
rp.button(panel4, action = showpos("c(100,100,90,20)"), title = "Button 2", pos = c(100,100,90,20))
rp.button(panel4, action = showpos("c(100,120,90,20"), title = "Button 3", pos = c(100,120,90,20))
rp.button(panel4, action = showpos("c(190,120,90,20"), title = "Button 4", pos = c(190,120,90,20))
rp.block(panel4)

rp.messagebox("You can also place objects on top of each other: objects are 'stacked' in order of creation, with the first created object at the 'bottom'.")

panel5 <- rp.control(title="Objects on top of each other", size=c(500,416))
rp.image(panel5, image.file, pos = c(0,0,500,416), id = "gulls.image", action = showpos("c(0,0,500,416)"))
rp.button(panel5, action = showpos("c(20,20,40,40"), title = "Button 1", pos = c(20,20,40,40))
rp.button(panel5, action = showpos("c(100,100,90,20)"), title = "Button 2", pos = c(100,100,90,20))
rp.button(panel5, action = showpos("c(145,112,90,20"), title = "Button 3", pos = c(145,112,90,20))
rp.button(panel5, action = showpos("c(190,124,90,20"), title = "Button 4", pos = c(190,124,90,20))
rp.block(panel5)

rp.messagebox("You can even mix the modes of operation, though that sometimes leads to unexpected effects.  Note, for example, how the 'size' specified is overridden (compare with the last example).")

panel6 <- rp.control(title="Mixed Methods", size=c(500,416))
rp.image(panel6, image.file, pos = "left", id = "gulls.image", action = showpos("'left'"))
rp.button(panel6, action = showpos("NULL"), title = "Button 1")
rp.button(panel6, action = showpos("c(100,100,90,20)"), title = "Button 2", pos = c(100,100,90,20))
rp.button(panel6, action = showpos("c(145,112,90,20"), title = "Button 3", pos = c(145,112,90,20))
rp.button(panel6, action = showpos("c(190,124,90,20"), title = "Button 4", pos = c(190,124,90,20))
rp.button(panel6, action = showpos("NULL"), title = "Button 5")
rp.button(panel6, action = showpos("NULL"), title = "Button 6")
rp.block(panel6)

rp.messagebox("rp.pos example done.")

