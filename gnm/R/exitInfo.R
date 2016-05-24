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

exitInfo <- function(object){
    conv <- object$converged
    if (conv)
        cat("Algorithm converged\n")
    else {
        cat("\nTolerance: ", object$tolerance, "\n")
        cat("\nAbsolute scores >= ",
            "tolerance * sqrt(tolerance + diag(information matrix)):\n\n")
        score <- abs(attr(conv, "score"))
        fail <- score >= attr(conv, "criterion")
        print(data.frame(abs.score = score,
                         criterion =  attr(conv, "criterion"))[fail,])
    }
}
