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

pogc <- function(oh, genotypeError)
{
    diag(oh) <- NA
    halfsib <- apply(oh, 1, function(x) names(which(x < genotypeError)))
    halfsib <- halfsib[unlist(lapply(halfsib, length)) > 2]
    if (length(halfsib) > 0)
    {
        repeated <- lapply(halfsib, length)
        pedigree <- data.frame(unlist(halfsib), rep(names(halfsib), repeated))
        colnames(pedigree) <- c("offspring", "parent")
        return(pedigree)
    }
    print("Could not find any groups!")
} 
