
## Library Loading
## Copyright 2010-2012 Jean-Pierre Rossi
## This file is part of the R-package `rich'.

.onLoad <- function(lib, pkg)  {
packageStartupMessage("This is rich ",
                           utils::packageDescription("rich", field="Version"),
                           appendLF = TRUE)}

