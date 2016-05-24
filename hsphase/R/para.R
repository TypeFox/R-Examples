# Copyright (C) 2014 Mohammad H. Ferdosi
#
# HSPhase is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# HSPhase program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http:#www.gnu.org/licenses/>.

para <- function(halfsibs, cpus = 1, option = "bmh", type = "SOCK", bmh_forwardVectorSize = 30, bmh_excludeFP = TRUE, bmh_nsap = 3, pmMethod = "constant")
{
    debut <- proc.time()
    if (!is.list(halfsibs)) 
        stop("halfsibs should be a list of matrices")
    bmhRnull <- function(x)
    {
        if (!is.null(x))
        {
            bmh(x, forwardVectorSize = bmh_forwardVectorSize, excludeFP = bmh_excludeFP, nsap = bmh_nsap)
        }
    }
    aioRnull <- function(x)
    {
        if (!is.null(x))
        {
            aio(x, bmh_forwardVectorSize = bmh_forwardVectorSize, bmh_excludeFP = bmh_excludeFP, bmh_nsap = bmh_nsap)
        }
    }
    rec <- function(x)
    {
        if (!is.null(x))
        {
            y <- pm(bmh(x, forwardVectorSize = bmh_forwardVectorSize, excludeFP = bmh_excludeFP, nsap = bmh_nsap), method = pmMethod)
            y <- apply(y, 2, sum)
            y[y > nrow(x)/2] <- 0
            # names(y)=colnames(x)[-1]
            y
        }
    }
    
    
    sfInit(parallel = TRUE, cpus = cpus, type = type)
    if (option == "bmh")
    {
        allhs <- sfLapply(halfsibs, bmhRnull)
    }
    if (option == "aio")
    {
        allhs <- sfLapply(halfsibs, aioRnull)
    }
    if (option == "rec")
    {
        allhs <- sfLapply(halfsibs, rec)
    }
    if (option == "pm")
    {
        allhs <- sfLapply(halfsibs, function(x) pm(bmh(x), pmMethod))
    }
    if (option == "ssp")
    {
        allhs <- sfLapply(halfsibs, function(x)
        {
            ssp(bmh(x, forwardVectorSize = bmh_forwardVectorSize, excludeFP = bmh_excludeFP, nsap = bmh_nsap), x)
        })
        sfRemoveAll()
        sfStop()
    }
    print(proc.time() - debut)
    allhs
} 
