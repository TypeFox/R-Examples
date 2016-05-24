draw <- function(panel) 
{
  plot(1:20, (1:20)^panel$h)
  panel
}
   
redraw <- function(panel) 
{
  rp.tkrreplot(panel, tkrp)
  panel
}

rpplot <- rp.control(title = "Demonstration of rp.tkrplot", h=1)
rp.tkrplot(rpplot, tkrp, draw)
rp.slider(rpplot, h, action=redraw, from=0.05, to=2.00, resolution=0.05)                                                                
