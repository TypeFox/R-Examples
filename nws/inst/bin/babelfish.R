#
# Copyright (c) 2005-2008, REvolution Computing, Inc.
#
# NetWorkSpaces is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
# USA
#

library(nws)

maxlen = 10000

bws = netWorkSpace('R babelfish')
repeat {
  s = paste(capture.output(nwsFetch(bws, 'food')), collapse = '\n')
  if (nchar(s) > maxlen) {
    s = paste(substr(s, 1, maxlen), '[WARNING: output truncated]')
  }
  nwsStore(bws, 'doof', s)
}
