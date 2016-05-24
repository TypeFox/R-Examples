.onAttach <- function(...) {

    mylib <- dirname(system.file(package = "WhatIf"))
    ver <- packageDescription("WhatIf", lib.loc = mylib)$Version
    builddate <- packageDescription("WhatIf", lib.loc = mylib)$Date

    # WhatIf Info - do not exceed 80char/line
 packageStartupMessage(paste("#######################################################\n",
"##\n",
"##  WhatIf (Version ", ver, ", built ", builddate, ")\n",
"##  Complete documentation available from http://gking.harvard.edu/whatif \n",
"##\n",
"#######################################################", sep=""))
}
