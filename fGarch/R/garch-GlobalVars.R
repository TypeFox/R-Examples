
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:               DESCRIPTION:
#  .setfGarchEnv           Set Global Vars
#  .garchArgsParser        Parses formula and data for garchFit
#  .garchOptimizerControl  Sets default values for Garch Optimizer
#  .garchFit               ... old Version, still in use by garchFit()
#  .garchNames             Slot names, @fit slot, parameters and controls
################################################################################


.setfGarchEnv(.llh = 1e99)
.setfGarchEnv(.garchDist = NA)
.setfGarchEnv(.params = NA)
.setfGarchEnv(.series = NA)
.setfGarchEnv(.trace = NA)  


################################################################################

