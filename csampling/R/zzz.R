.onAttach <- function(libname, pkgname) 
{
  version <- as.character("1.2-2 (2014-03-31)")
  psm <- paste("\n   Package \"csampling\"", version, "\n",
    "   Copyright (C) 2000-2014 A. R. Brazzale\n\n",
    "This is free software, and you are welcome to redistribute\n",
    "it and/or modify it under the terms of the GNU General\n",
    "Public License published by the Free Software Foundation.\n",
    "Package \"csampling\" comes with ABSOLUTELY NO WARRANTY.\n\n",
    "type `help(package=\"csampling\")' for summary information\n")
  packageStartupMessage(psm)
  invisible()
}
