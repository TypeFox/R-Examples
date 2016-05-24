panel <- rp.control(x=1)
callback <- function(panel)
{
  rp.text.change(panel, "t2", panel$x)
  panel$x = panel$x+1
  panel
}
rp.text(panel, "This is a test", name="t1")
rp.text(panel ,"And so is this", font="Arial", foreground="white", background="navy", action=callback, name="t2")
rp.text(panel,"Here is some more text, this time across several lines.\nHere is some more text, this time across several lines.\nHere is some more text, this time across several lines.", name="t3")
