##########################################################################
#                                                                        #
#  SPRINT: Simple Parallel R INTerface                                   #
#  Copyright © 2008-2011 The University of Edinburgh                     #
#                                                                        #
#  This program is free software: you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  any later version.                                                    #
#                                                                        #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
#  GNU General Public License for more details.                          #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with this program. If not, see <http://www.gnu.or/licenses/>.   #
#                                                                        #
##########################################################################

init.rng <- function(seed)
{
  if ( missing(seed) | is.null(seed) ) {
    seed <- .Call("sprint_create_seed")
  }
  ## Fix up for max value allowed by rlecuyer
  seed <- (seed %% 4294944442) + 1

  invisible(.Call("init_rng", seed))
}

init.rng.worker <- function(rank, size, seed)
{
  name <- rank + 1

  .lec.SetPackageSeed(rep(seed, 6))

  .lec.CreateStream(1:size)

  kind.old <- .lec.CurrentStream(name)

  return(kind.old)
}

reset.rng <- function()
{
  invisible(.Call("reset_rng"))
}

reset.rng.worker <- function(kind.old, size)
{
  .lec.CurrentStreamEnd(kind.old)
  .lec.DeleteStream(1:size)
}
