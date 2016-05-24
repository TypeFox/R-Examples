.onAttach <- function(...) {
  if (!interactive()){return()}

  the_startup_message<-paste(sep="",
  "This is HotDeckImputation version ",utils::packageVersion("HotDeckImputation"),".\n",
  "Academic users, please be sure to use: citation(\"HotDeckImputation\").\n",
  "Feel free to contact the Maintainer about licensing, feature requests, etc.\n",
  "Use suppressPackageStartupMessages to eliminate package startup messages.")
  
  packageStartupMessage(the_startup_message)
}
