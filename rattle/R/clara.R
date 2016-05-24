# Gnome R Data Miner: GNOME interface to R for Data Mining
#
# Time-stamp: <2011-06-23 21:17:19 Graham Williams>
#
# Implement biclust functionality.
#
# Copyright (c) 2010 Togaware Pty Ltd
#
# This files is part of Rattle.
#
# Rattle is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# Rattle is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rattle. If not, see <http://www.gnu.org/licenses/>.

########################################################################
# ToDo 100121
#
# Execute.
# Graphical display of output.
# Allow choice of methods.

########################################################################
# Callbacks

# When a radio button is selected, display the appropriate tab page.

on_clara_radiobutton_toggled <- function(button)
{
  if (button$getActive())
    crv$CLUSTER$setCurrentPage(crv$CLUSTER.CLARA.TAB)
  setStatusBar()
}

