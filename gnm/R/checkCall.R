#  Copyright (C) 2006 Heather Turner
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

checkCall <- function(){
    badCall <- lapply(sys.calls(), "[[", 1) %in% c("model.frame.default",
                                                   "model.matrix.default")
    if (any(badCall))
        stop(paste(sys.call(-1)[[1]], "terms are only valid in gnm models."))
}
