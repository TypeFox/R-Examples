###########################################################################
# Copyright 2009 Nobody                                                   #
#                                                                         #
# This file is part of hergm.                                             #
#                                                                         # 
#    hergm is free software: you can redistribute it and/or modify        #
#    it under the terms of the GNU General Public License as published by #
#    the Free Software Foundation, either version 3 of the License, or    #
#    (at your option) any later version.                                  #
#                                                                         # 
#    hergm is distributed in the hope that it will be useful,             #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
#    GNU General Public License for more details.                         #
#                                                                         #
#    You should have received a copy of the GNU General Public License    #
#    along with hergm.  If not, see <http://www.gnu.org/licenses/>.       #
#                                                                         # 
###########################################################################

.onAttach <- function(lib, pkg)
{
  info <- packageDescription("hergm")
  packageStartupMessage(
  paste('\nhergm: version ', info$Version, ', created on ', info$Date, '\n',
          "Copyright (c) 2009, Michael Schweinberger, Rice University\n", 
          "                    with contributions from:\n",
          "                    Mark S. Handcock, University of California - Los Angeles\n",
          "                    Pamela Luna, Rice University\n",
          "To start hergm: enter help(package=\"hergm\")\n",
          'For license and citation information type citation("hergm").\n', sep="")
 )
}

