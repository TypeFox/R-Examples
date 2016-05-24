## Open html help file

openHelpFile <-
function(filename){
  print(help(filename, package = "DALY", help_type = "html"))
}