callback <- function(panel)
{
  print(panel$option)
}

panel <- rp.control()
rp.combo(panel, option, "Pick an option:", c("Option1","Option2","Other options"), action=callback)

