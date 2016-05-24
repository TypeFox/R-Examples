##
## Interactive graphical elicitation of split-normal parameters for the UK future net migration.
##
library("fanplot")

##
##install packages if not already done so (uncomment)
##
if (!(require("shiny"))) {
  install.packages("shiny")
  library("shiny")
}

##
##run shiny app
##
runApp(system.file("netelicit", package = "fanplot"))