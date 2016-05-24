callback <- function(panel)
{
  print(panel$option)
}

panel1 <- rp.control()
rp.textentry(panel1, option, labels="Your name:", initval="-", action=callback, width=40)

panel2 <- rp_window()
rp.textentry(panel2, option, labels=c("Your height:", "Your weight:"), initval=c("H", "W"), action=callback, width=20)

