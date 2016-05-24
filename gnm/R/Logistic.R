#  Copyright (C) 2007 Heather Turner
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

Logistic <- function(x, inst = NULL){
    list(predictors = list(Asym = 1, xmid = 1, scal = 1),
         variables = list(substitute(x)),
         term = function(predLabels, varLabels) {
             paste(predLabels[1], "/(1 + exp((", predLabels[2], "-",
                   varLabels[1], ")/", predLabels[3], "))")
         },
         start = function(theta){
             c(NA, mean(x), sd(x))
         }
         )
}
class(Logistic) <- "nonlin"
