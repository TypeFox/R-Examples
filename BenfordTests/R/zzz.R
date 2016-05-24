.onAttach <- function(...) {
  if (!interactive()){return()}

  the_startup_message<-paste(sep="",
  "This is BenfordTests version ",utils::packageVersion("BenfordTests"),".\n",
  "Academic users, please be sure to use: citation(\"BenfordTests\").\n",
  "Feel free to contact the Maintainer about liscensing, feature requests, etc.\n",
  "Use suppressPackageStartupMessages to eliminate package startup messages.")
  
  packageStartupMessage(the_startup_message)
}
