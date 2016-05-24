
### Copyright (C) 2001-2006  Deepayan Sarkar <Deepayan.Sarkar@R-project.org>
### Copyright (C) 2001-2005  Saikat DebRoy <saikat@stat.wisc.edu>
###
### This file is part of the lattice package for R.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
### MA 02110-1301, USA




.cupdate <- function(index, maxim)
{
	
	## This unexported function is used to handle arbitrary number of
	## conditioning variables : every time it is called, it increments
	## the "current" level of the conditioning variables suitably,
	## i.e., it tries to increment the level of the 1st conditining
	## variable (the one which varies fastest along panel order) and
	## if it happens to be at its maximum (last) value, it sets it to
	## the first value AND increments the "current" level of the 2nd
	## (next) conditioning variable recursively.
	
	if(length(index)!=length(maxim)||length(maxim)<=0)
		stop("Inappropriate arguments")
	index[1] <- index[1] + 1
	if (index[1] > maxim[1] && length(maxim) > 1)
		c(1, .cupdate(index[-1], maxim[-1]))
	else index
}



.check.layout <-
    function(layout, cond.max.level, skip = FALSE)
{
    if (all(skip)) stop("skip cannot be all TRUE")
    number.of.cond <- length(cond.max.level)
    nplots <- prod(cond.max.level)

    if (!is.numeric(layout))
    {
        layout <- c(0,1,1)
        if (number.of.cond == 1) layout[2] <- nplots
        else
        {
            layout[1] <- cond.max.level[1]
            layout[2] <- cond.max.level[2]
        }
        skip <- rep(skip, length.out = max(layout[1] * layout[2], layout[2]))
        plots.per.page <- length(skip) - length(skip[skip])
        layout[3] <- ceiling(nplots/plots.per.page) # + 1
    }
    else if (length(layout) == 1)
        stop("layout must have at least 2 elements")
    else if (length(layout) == 2)
    {
        if (all(is.na(layout)))
            stop("inadmissible value of layout")
        else if (all(layout < 1))
            stop("at least one element of layout must be positive")
        else if (isTRUE(layout[2] == 0))
            stop("inadmissible value of layout")

        if (is.na(layout[1]))
            layout[1] <- ceiling(nplots / layout[2])
        if (is.na(layout[2]))
            layout[2] <- ceiling(nplots / layout[1])

        skip <- rep(skip, length.out = max(layout[1] * layout[2], layout[2]))
        plots.per.page <- length(skip) - length(skip[skip])
        layout[3] <- ceiling(nplots / plots.per.page) # + 1
    }
    else if (length(layout)==3)
    {
        if(layout[1] < 0 || layout[2] < 1 || layout[3] < 1)
            stop("invalid value for layout")
    }
    layout
}


.compute.packet <-
    function(cond, levels)
{
    id <- !(do.call("pmax", lapply(cond, is.na)))
    stopifnot(any(id))
    for (i in seq_along(cond))
    {
        var <- cond[[i]]
        id <-  id & (as.numeric(var) == levels[i])
    	## MARCO: Removed the possibility of numerical conditioning variables
	}
    id
}

