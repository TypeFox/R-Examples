#  R package rjags file R/windows/progress.bar.R
#  Copyright (C) 2009 Martyn Plummer
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License version
#  2 as published by the Free Software Foundation.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

updatePB <- function(start.iter, end.iter, adapting)
{
    winProgressBar(title=ifelse(adapting, "Adapting","Updating"),
                  label="Iteration 0", min = start.iter, max=end.iter,
                  initial=start.iter)
}

setPB <- function(pb, iter)
{
    setWinProgressBar(pb, iter, label=paste("Iteration",iter))
}
