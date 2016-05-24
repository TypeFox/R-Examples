.onAttach <- function(libname, pkgname) 
{
  version <- as.character("1.2-3 (2014-06-27)")
  psm <- paste("\n   Package \"cond\"", version, "\n", 
    "   Copyright (C) 2000-2014 A. R. Brazzale\n\n",
    "This is free software, and you are welcome to redistribute\n",
    "it and/or modify it under the terms of the GNU General\n",
    "Public License published by the Free Software Foundation.\n",
    "Package \"cond\" comes with ABSOLUTELY NO WARRANTY.\n\n",
    "type `help(package=\"cond\")' for summary information\n")
  packageStartupMessage(psm)
  invisible()
}
