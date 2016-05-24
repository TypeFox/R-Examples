## Copyright (C) Yusman Kamaleri

## This file is part of aimPlot.

## require(aimPlot)
## detach("package:aimPlot")

.onAttach <- function(lib, pkg,...){
   packageStartupMessage(installrWelcomeMessage())
}


installrWelcomeMessage <- function(){
   
   paste("\n",     
         "Welcome to aimPlot version ", utils::packageDescription("aimPlot")$Version, "\n",
         "\n",
         # "Type ?aimPlot to access the overall documentation and\n",
         "More information is available on the Github:\n",
         "https://github.com/ybkamaleri/aimPlot/\n",
         "\n",               
         "Contact: <ybkamaleri@gmail.com>\n",
         "Suggestions and bug-reports can be submitted at: https://github.com/ybkamaleri/aimPlot/issues\n",
         "\n",
         "\tTo suppress this message use:\n",
         "\tsuppressPackageStartupMessages(library(aimPlot))\n",  
         sep="")
}
