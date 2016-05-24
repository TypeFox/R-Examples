#example function for "pos" documentation.
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

rp.pos <- function(layout = "default") {

if (layout == "default")
{
  rp.messagebox("This is a demonstration of pos using its default (unset) value. Resize the windows to see the behaviour.")
  # Create an rpanel and add some buttons to it: 'default' mode

  panel1 <- rp.control(title='Default mode: no pos specified')
  rp.button(panel1, action = showpos('NULL'), title = "Button 1")
  rp.button(panel1, action = showpos('NULL'), title = "Button 2")
  rp.button(panel1, action = showpos('NULL'), title = "Button 3")
}

if (layout == "pack")
{
  rp.messagebox("This is a demonstration of pos using its 'pack' mode. Resize the windows to see the behaviour.")

  #create an rpanel and put buttons in using "pack" options.
  panel2 <- rp.control(title="Pack mode: pos=left/right/top/bottom")
  rp.button(panel2, action = showpos("'left'"), title = "Button 1", pos = "left")
  rp.button(panel2, action = showpos("'left'"), title = "Button 2", pos = "left")
  rp.button(panel2, action = showpos("'right'"), title = "Button 3", pos = "right")
  rp.button(panel2, action = showpos("'top'"), title = "Button 4", pos = "top")
  rp.button(panel2, action = showpos("'bottom'"), title = "Button 5", pos = "bottom")
  rp.button(panel2, action = showpos("'right'"), title = "Button 6", pos = "right")
  rp.button(panel2, action = showpos("'bottom'"), title = "Button 7", pos = "bottom")
  # rp.block(panel2)

  rp.messagebox("A more typical use of this is to have an interactive image or a Tkrplot to the side of some buttons...")

  panel3 <- rp.control(title="Graphic and buttons, using 'pack'")
  image.file <- file.path(system.file(package = "rpanel"), "images", "gulllmks.gif")
  rp.image(panel3, image.file, pos = "left", name = "gulls.image", action = showpos("'left'"))
  rp.button(panel3, action = showpos("'top'"), title = "Button 1", pos = "top")
  rp.button(panel3, action = showpos("'top'"), title = "Button 2", pos = "top")
  rp.button(panel3, action = showpos("'top'"), title = "Button 3", pos = "top")
  rp.button(panel3, action = showpos("'top'"), title = "Button 4", pos = "top")
  rp.button(panel3, action = showpos("'bottom'"), title = "Button 5", pos = "bottom")
  # rp.block(panel3)
}

if (layout == "place")
{
  rp.messagebox("Specifying pixel positions allows a huge amount of control, but takes more time to do and does not resize dynamically")

  panel4 <- rp.control(title="'place' demo",size=c(300,150))
  rp.button(panel4, action = showpos("c(20,20,40,40"), title = "Button 1", pos = c(20,20,40,40))
  rp.button(panel4, action = showpos("c(100,100,90,20)"), title = "Button 2", pos = c(100,100,90,20))
  rp.button(panel4, action = showpos("c(100,120,90,20"), title = "Button 3", pos = c(100,120,90,20))
  rp.button(panel4, action = showpos("c(190,120,90,20"), title = "Button 4", pos = c(190,120,90,20))
  # rp.block(panel4)

  rp.messagebox("You can also place objects on top of each other: objects are 'stacked' in order of creation, with the first created object at the 'bottom'.")

  panel5 <- rp.control(title="Objects on top of each other", size=c(500,416))
  image.file <- file.path(system.file(package = "rpanel"), "images", "gulllmks.gif")
  rp.image(panel5, image.file, pos = c(0,0,500,416), name = "gulls.image", action = showpos("c(0,0,500,416)"))
  rp.button(panel5, action = showpos("c(20,20,40,40"), title = "Button 1", pos = c(20,20,40,40))
  rp.button(panel5, action = showpos("c(100,100,90,20)"), title = "Button 2", pos = c(100,100,90,20))
  rp.button(panel5, action = showpos("c(145,112,90,20"), title = "Button 3", pos = c(145,112,90,20))
  rp.button(panel5, action = showpos("c(190,124,90,20"), title = "Button 4", pos = c(190,124,90,20))
  # rp.block(panel5)
}

if ((layout == "place") || (layout == "pack"))
{
  rp.messagebox("You can even mix the 'pack' and 'place' modes of operation, though that sometimes leads to unexpected effects.  Note, for example, how the 'size' specified is overridden (compare with the last example).")

  panel6 <- rp.control(title="Mixed Methods", size=c(500,416))
  image.file <- file.path(system.file(package = "rpanel"), "images", "gulllmks.gif")
  rp.image(panel6, image.file, pos = "left", name = "gulls.image", action = showpos("'left'"))
  rp.button(panel6, action = showpos("NULL"), title = "Button 1")
  rp.button(panel6, action = showpos("c(100,100,90,20)"), title = "Button 2", pos = c(100,100,90,20))
  rp.button(panel6, action = showpos("c(145,112,90,20"), title = "Button 3", pos = c(145,112,90,20))
  rp.button(panel6, action = showpos("c(190,124,90,20"), title = "Button 4", pos = c(190,124,90,20))
  rp.button(panel6, action = showpos("NULL"), title = "Button 5")
  rp.button(panel6, action = showpos("NULL"), title = "Button 6")
  # rp.block(panel6)
}

if (layout == "grid")
{
  rp.messagebox("This is an unlikely demonstration of grid. rp.geosim, rp.mururoa and rp.rosyth are good examples for its use.")

  pnl <- rp.control()
  rp.grid(pnl, "g1", pos=list(row=0, column=0), background="red")
  rp.grid(pnl, "g2", pos=list(row=1, column=0), background="navy")
  rp.grid(pnl, "g3", pos=list(row=2, column=0), background="green")
  rp.grid(pnl, "g4", pos=list(row=0, column=1, rowspan=4), background="yellow")
  rp.button(pnl, function(panel) { panel}, "Press me 1", pos= list(row=0,column=0,grid="g1"))
  rp.button(pnl, function(panel) { panel}, "Press me 2", pos= list(row=1,column=1,grid="g1"))
  rp.button(pnl, function(panel) { panel}, "Press me 3", pos= list(row=2,column=2,grid="g1"))
  rp.button(pnl, function(panel) { panel}, "Press me 4", pos= list(row=0,column=0,grid="g2"))
  rp.button(pnl, function(panel) { panel}, "Press me 5", pos= list(row=1,column=1,grid="g2"))
  rp.button(pnl, function(panel) { panel}, "Press me 6", pos= list(row=2,column=2,grid="g2"))
  rp.button(pnl, function(panel) { panel}, "Press me 7", pos= list(row=0,column=0,grid="g3"))
  rp.button(pnl, function(panel) { panel}, "Press me 8", pos= list(row=1,column=1,grid="g3"))
  rp.button(pnl, function(panel) { panel}, "Press me 9", pos= list(row=2,column=2,grid="g3"))
  rp.button(pnl, function(panel) { panel}, "Press me 10", pos= list(row=0,column=0,grid="g4"))
  rp.button(pnl, function(panel) { panel}, "Press me 11", pos= list(row=1,column=1,grid="g4"))
  rp.button(pnl, function(panel) { panel}, "Press me 12", pos= list(row=2,column=2,grid="g4"))
  rp.button(pnl, function(panel) { panel}, "Press me 10", pos= list(row=3,column=0,grid="g4"))
  rp.button(pnl, function(panel) { panel}, "Press me 11", pos= list(row=4,column=1,grid="g4"))
  rp.button(pnl, function(panel) { panel}, "Press me 12", pos= list(row=5,column=2,grid="g4"))
}

rp.messagebox("Demonstration complete. Please close all the rpanel windows.")
invisible()

}

